package ahrd.model;

public class Fscore {
	private Double score = 0.0;
	private Double precision = 0.0; // positive predictive value (PPV) = TP / prediction 
	private Double recall = 0.0; // true positive rate (TPR) = TP / groundTruth
	private Double remainingUncertainty = 0.0; // false discovery rate (FDR) = FP / prediction = 1 - PPV
	private Double missingInformation = 0.0; // false negative rate (FNR) = FN / groundTruth = 1 - TPR

	public Fscore() {
		super();
	}
	
	public Fscore(double s, double p, double r) {
		this.setScore(s);
		this.setPrecision(p);
		this.setRecall(r);
	}

	public Double getScore() {
		return score;
	}

	public void setScore(Double score) {
		this.score = score;
	}

	public Double getPrecision() {
		return precision;
	}

	public void setPrecision(Double precision) {
		this.precision = precision;
	}

	public Double getRecall() {
		return recall;
	}

	public void setRecall(Double recall) {
		this.recall = recall;
	}

	public Double getRemainingUncertainty() {
		return remainingUncertainty;
	}

	public void setRemainingUncertainty(Double remainingUncertainty) {
		this.remainingUncertainty = remainingUncertainty;
	}

	public Double getMissingInformation() {
		return missingInformation;
	}

	public void setMissingInformation(Double missingInformation) {
		this.missingInformation = missingInformation;
	}

}
