package ahrd.model;

public class InfoContentFscore extends Fscore{
	private Double commonInfoContentPrediction = 0.0;
	private Double infoContentGroundTruth = 0.0;
	private Double commonInfoContentGroundTruth = 0.0;
	private Double infoContentPrediction = 0.0;
	
	public InfoContentFscore() {
		super();
	}

	public Double getCommonInfoContentPrediction() {
		return commonInfoContentPrediction;
	}

	public void setCommonInfoContentPrediction(Double commonInfoContentPrediction) {
		this.commonInfoContentPrediction = commonInfoContentPrediction;
	}

	public Double getInfoContentGroundTruth() {
		return infoContentGroundTruth;
	}

	public void setInfoContentGroundTruth(Double infoContentGroundTruth) {
		this.infoContentGroundTruth = infoContentGroundTruth;
	}

	public Double getCommonInfoContentGroundTruth() {
		return commonInfoContentGroundTruth;
	}

	public void setCommonInfoContentGroundTruth(Double commonInfoContentGroundTruth) {
		this.commonInfoContentGroundTruth = commonInfoContentGroundTruth;
	}

	public Double getInfoContentPrediction() {
		return infoContentPrediction;
	}

	public void setInfoContentPrediction(Double infoContentPrediction) {
		this.infoContentPrediction = infoContentPrediction;
	}
}
