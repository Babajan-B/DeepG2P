from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile

from fastapi import FastAPI, File, HTTPException, UploadFile
from pydantic import BaseModel

from deepg2p_celiac.api.prediction import load_attribution, predict_graph_smoothed, predict_snp_panel

APP_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PANEL = APP_ROOT / "configs" / "celiac_snp_panel.csv"
DEFAULT_RUN = APP_ROOT / "data" / "processed" / "run"
DEFAULT_MODEL = DEFAULT_RUN / "model" / "best_model.joblib"
DEFAULT_GRAPH_MODEL = DEFAULT_RUN / "graph_model" / "best_graph_model.joblib"
DEFAULT_ADJACENCY = DEFAULT_RUN / "graph" / "gene_adjacency.csv"
DEFAULT_TOP_SNPS = DEFAULT_RUN / "model" / "top_snp_features.csv"
DEFAULT_TOP_GENES = DEFAULT_RUN / "graph_model" / "top_gene_features.csv"

app = FastAPI(
    title="DeepG2P-Celiac",
    description="Research-only celiac disease SNP-panel prediction API.",
    version="0.1.0",
)


class PredictionResponse(BaseModel):
    sample_id: str
    model_name: str
    research_risk_score: float
    predicted_label: str
    top_features: list[dict] = []
    warning: str


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/artifacts")
def artifacts() -> dict[str, bool]:
    return {
        "panel": DEFAULT_PANEL.exists(),
        "snp_model": DEFAULT_MODEL.exists(),
        "graph_model": DEFAULT_GRAPH_MODEL.exists(),
        "graph_adjacency": DEFAULT_ADJACENCY.exists(),
        "top_snps": DEFAULT_TOP_SNPS.exists(),
        "top_genes": DEFAULT_TOP_GENES.exists(),
    }


@app.post("/predict", response_model=PredictionResponse)
async def predict(file: UploadFile = File(...)) -> PredictionResponse:
    if not DEFAULT_MODEL.exists():
        raise HTTPException(status_code=503, detail="No trained SNP-panel model found.")

    suffix = Path(file.filename or "sample.txt").suffix or ".txt"
    with NamedTemporaryFile("wb", suffix=suffix, delete=True) as tmp:
        tmp.write(await file.read())
        tmp.flush()
        result = predict_snp_panel(Path(tmp.name), DEFAULT_PANEL, DEFAULT_MODEL)

    return PredictionResponse(
        **result,
        top_features=load_attribution(DEFAULT_TOP_SNPS),
        warning="Research-only estimate. Not a medical diagnosis.",
    )


@app.post("/predict/graph", response_model=PredictionResponse)
async def predict_graph(file: UploadFile = File(...)) -> PredictionResponse:
    missing = [
        name
        for name, path in {
            "graph model": DEFAULT_GRAPH_MODEL,
            "graph adjacency": DEFAULT_ADJACENCY,
        }.items()
        if not path.exists()
    ]
    if missing:
        raise HTTPException(status_code=503, detail=f"Missing artifacts: {', '.join(missing)}")

    suffix = Path(file.filename or "sample.txt").suffix or ".txt"
    with NamedTemporaryFile("wb", suffix=suffix, delete=True) as tmp:
        tmp.write(await file.read())
        tmp.flush()
        result = predict_graph_smoothed(
            Path(tmp.name),
            DEFAULT_PANEL,
            DEFAULT_GRAPH_MODEL,
            DEFAULT_ADJACENCY,
        )

    return PredictionResponse(
        **result,
        top_features=load_attribution(DEFAULT_TOP_GENES),
        warning="Research-only estimate. Not a medical diagnosis.",
    )
