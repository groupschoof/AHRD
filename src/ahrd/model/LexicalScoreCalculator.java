package ahrd.model;


public class LexicalScoreCalculator {

	private Protein protein;

	public LexicalScoreCalculator(Protein protein) {
		setProtein(protein);
	}

	public double sumTokenScoresDividebByHighScore(BlastResult br) {
		TokenScoreCalculator tsc = getProtein().getTokenScoreCalculator();
		return tsc.sumOfAllTokenScores(br) / tsc.getTokenHighScore();
	}

	public double lexicalScore(BlastResult br) {
		return ((sumTokenScoresDividebByHighScore(br) / correctionFactor(br)));
	}

	/**
	 * For argument BlastResult: Number of br's tokens / Number of br's
	 * informative tokens
	 */
	public double correctionFactor(BlastResult br) {
		TokenScoreCalculator tsc = getProtein().getTokenScoreCalculator();
		double noInformativeTokens = 0.0;
		for (String token : br.getTokens()) {
			if (tsc.isInformativeToken(token))
				noInformativeTokens += 1.0;
		}
		return ((double) br.getTokens().size()) / noInformativeTokens;
	}

	public Protein getProtein() {
		return protein;
	}

	public void setProtein(Protein protein) {
		this.protein = protein;
	}

}
