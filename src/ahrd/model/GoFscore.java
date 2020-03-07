package ahrd.model;

public class GoFscore extends Fscore{
	private Fscore allFscore;
	private Fscore bpoFscore;
	private Fscore mfoFscore;
	private Fscore ccoFscore;

	public GoFscore() {
		super();
	}

	@Override
	public Double getPrecision() {
		return allFscore.getPrecision();
	}

	@Override
	public void setPrecision(Double precision) {
		allFscore.setPrecision(precision);
		allFscore.calcScore();
	}

	@Override
	public Double getRecall() {
		return allFscore.getRecall();
	}

	@Override
	public void setRecall(Double recall) {
		allFscore.setRecall(recall);
		allFscore.calcScore();
	}

	@Override
	public Double getScore() {
		return allFscore.getScore();
	}

	public Fscore getAllFscore() {
		return allFscore;
	}

	public void setAllFscore(Fscore allFscore) {
		this.allFscore = allFscore;
	}

	public Fscore getBpoFscore() {
		return bpoFscore;
	}

	public void setBpoFscore(Fscore bpoFscore) {
		this.bpoFscore = bpoFscore;
	}

	public Fscore getMfoFscore() {
		return mfoFscore;
	}

	public void setMfoFscore(Fscore mfoFscore) {
		this.mfoFscore = mfoFscore;
	}

	public Fscore getCcoFscore() {
		return ccoFscore;
	}

	public void setCcoFscore(Fscore ccoFscore) {
		this.ccoFscore = ccoFscore;
	}

}
