import scanpy as sc
import scvi
import ray
from ray import tune
from ray.tune.schedulers import ASHAScheduler
import torch
import sys, os
import matplotlib.pyplot as plt
# Ray 초기화 - 최대 16개의 CPU 스레드 사용하도록 설정
ray.init(num_cpus=56)

# 데이터 읽기
adata = sc.read(sys.argv[1])
outputdir = sys.argv[2]
os.makedirs(outputdir, exist_ok=True)
if os.path.exists(outputdir + '/the_model/'): print('"'+outputdir + '/the_model/" Already exists'); exit()
# 유전자 필터링
sc.pp.filter_genes(adata, min_cells=100)

# 모델 설정
model_cls = scvi.model.SCVI

# 하이퍼파라미터 검색 공간 설정
search_space = {
	"n_hidden": tune.choice([92, 128, 192, 256]),
	"n_latent": tune.choice([10, 20, 30, 40, 50, 60]),
	"n_layers": tune.choice([1, 2, 3]),
	"lr": tune.loguniform(1e-4, 1e-2),
	"gene_likelihood": tune.choice(["nb", "zinb"])}

# 훈련 함수 정의
def train_scvi(config, adata):
	scvi.model.SCVI.setup_anndata(adata, categorical_covariate_keys=['Sample', 'DX'], continuous_covariate_keys=['pct_counts_mt', 'pct_counts_ribo'])
	model = scvi.model.SCVI(adata, n_hidden=config["n_hidden"], n_latent=config["n_latent"], n_layers=config["n_layers"], gene_likelihood=config["gene_likelihood"])
	# 학습 속도 설정
	model.optimizer_cls = lambda x: torch.optim.Adam(x, lr=config["lr"])
	model.train(max_epochs=20)
	# validation_loss를 단일 값으로 반환
	validation_loss = model.history["elbo_train"].iloc[-1].item()
	return {"validation_loss": validation_loss}
# 튜닝 실행
scheduler = ASHAScheduler(metric="validation_loss", mode="min")
analysis = tune.run(
	tune.with_parameters(train_scvi, adata=adata),
	resources_per_trial={"cpu": 14, "gpu": 0},
	config=search_space,num_samples=100,
	scheduler=scheduler
)

# 결과 확인
results_df = analysis.results_df
# Ray 종료
ray.shutdown()

results_df.to_csv(outputdir + "tune_results.csv", index=False)

best_result = results_df.loc[results_df['validation_loss'].idxmin()]
print(best_result)

n_hidden, n_latent, n_layers, lr, gene_likelihood = best_result['config/n_hidden'], best_result['config/n_latent'], best_result['config/n_layers'], best_result['config/lr'], best_result['config/gene_likelihood']

# scVI Integration
scvi.model.SCVI.setup_anndata(adata, categorical_covariate_keys = ['Sample', 'DX'], continuous_covariate_keys=['pct_counts_mt', 'pct_counts_ribo'])



model = scvi.model.SCVI(adata, n_hidden = n_hidden, n_latent = n_latent, n_layers = n_layers, gene_likelihood = gene_likelihood)
kwargs = {'lr':lr}

model.train(max_epochs = 200, early_stopping = True, plan_kwargs = kwargs)

model.save(outputdir + '/the_model')

y = model.history['reconstruction_loss_validation']['reconstruction_loss_validation'].min()
adata.write_h5ad(outputdir + '/temp.h5ad')

plt.plot(model.history['reconstruction_loss_train']['reconstruction_loss_train'], label='train')
plt.plot(model.history['reconstruction_loss_validation']['reconstruction_loss_validation'], label='validation')

plt.axhline(y, c = 'k')

plt.legend()
plt.show()
plt.savefig(outputdir + '/train_validation.png', bbox_inches='tight')
