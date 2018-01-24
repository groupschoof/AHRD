package ahrd.model;

import static ahrd.controller.Settings.getSettings;

public class Fscore {
	private Double score = 0.0;
	private Double precision = 0.0; // positive predictive value (PPV) = TP / prediction 
	private Double recall = 0.0; // true positive rate (TPR) = TP / groundTruth
	private Double remainingUncertainty = 1.0; // false discovery rate (FDR) = FP / prediction = 1 - PPV
	private Double missingInformation = 1.0; // false negative rate (FNR) = FN / groundTruth = 1 - TPR

	public Fscore() {
		super();
	}
	
	public Fscore(double s, double p, double r) {
		this.score = s;
		this.precision = p;
		this.recall = r;
	}

	public Double getScore() {
		return score;
	}

	// F-Beta-Measure is the harmonic mean of precision and recall weighted by param beta:
	private void calcScore() {
		if (this.precision > 0.0 && this.recall > 0.0) {
			Double bSqr = getSettings().getFMeasureBetaParameter() * getSettings().getFMeasureBetaParameter();
			this.score = (1 + bSqr) * (this.precision * this.recall) / (bSqr * this.precision + this.recall);
		} else {
			this.score = 0.0;
		}
	}

	public Double getPrecision() {
		return precision;
	}

	public void setPrecision(Double precision) {
		this.precision = precision;
		this.remainingUncertainty = 1-precision;
		this.calcScore();
	}

	public Double getRecall() {
		return recall;
	}

	public void setRecall(Double recall) {
		this.recall = recall;
		this.missingInformation = 1 - recall;
		this.calcScore();
	}

	public Double getRemainingUncertainty() {
		return remainingUncertainty;
	}

	public Double getMissingInformation() {
		return missingInformation;
	}

}
