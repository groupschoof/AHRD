package ahrd.model;

public class Fscore {
	private Double score = 0.0;
	private Double precision = 0.0;
	private Double recall = 0.0;

	public Fscore() {
		super();
	}
	
	public Fscore(double s, double p, double r) {
		this.setScore(s);
		this.setPrecision(p);
		this.setScore(r);
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

}
