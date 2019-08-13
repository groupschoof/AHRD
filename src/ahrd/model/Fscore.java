package ahrd.model;

import static ahrd.controller.Settings.getSettings;

public class Fscore {
	private Double score = 0.0;
	private Double precision = 0.0; // positive predictive value (PPV) = TP / prediction 
	private Double recall = 0.0;  // true positive rate (TPR) = TP / groundTruth

	public Fscore() {
		super();
	}
	
	/**
	 * Should only be used in unit test classes
	 * 
	 * @param score
	 * @param precision
	 * @param recall
	 */
	public Fscore(double s, double p, double r) {
		this.score = s;
		this.precision = p;
		this.recall = r;
	}

	public Double getScore() {
		return score;
	}
	
	// F-Beta-Measure is the harmonic mean of precision and recall weighted by the beta parameter:
	// Arithmetically the fscore for a precision of 0 and a recall of 0 is undefined and should be stored as NaN (1 * 0 * 0 / (1 * 0 + 0) = 0/0 = NaN).
	// When no predictions are made the precision is NaN. When no ground Truth exists the recall is NaN. Both cases result in an NaN-fscore as well.
	// To differentiate these cases from a genuine prediction that has no overlap with a proper ground truth (p=0 and r=0) the fscore is set to 0 instead of NaN. 
	private void calcScore() {
		if (this.precision == 0.0 && this.recall == 0.0) {
			this.score = 0.0;
		} else {
			Double bSqr = getSettings().getFMeasureBetaParameter() * getSettings().getFMeasureBetaParameter();
			this.score = (1 + bSqr) * (this.precision * this.recall) / (bSqr * this.precision + this.recall);
		}
	}

	public Double getPrecision() {
		return precision;
	}

	public void setPrecision(Double precision) {
		this.precision = precision;
		this.calcScore();
	}

	public Double getRecall() {
		return recall;
	}

	public void setRecall(Double recall) {
		this.recall = recall;
		this.calcScore();
	}

}
