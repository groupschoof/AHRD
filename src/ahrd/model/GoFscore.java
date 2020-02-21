package ahrd.model;

public class GoFscore extends Fscore{
	private Fscore bpoFscore;
	private Fscore mfoFscore;
	private Fscore ccoFscore;

	public GoFscore() {
		super();
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
