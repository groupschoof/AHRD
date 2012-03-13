package ahrd.model;

import java.util.regex.Pattern;

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
		return ((sumTokenScoresDividebByHighScore(br) / correctionFactor(br) + geneOntologyScore(br)));
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
		return (new Double(br.getTokens().size()).doubleValue() / noInformativeTokens);
	}

	/**
	 * For each blastResult's token find GeneOntologyResults, where the token
	 * appears in the GO's name. Return the sum of those GeneOntologyResults'
	 * probabilities as geneOntologyScore.
	 */
	public double geneOntologyScore(BlastResult blastResult) {
		double goScore = 0.0;
		for (String token : blastResult.getTokens()) {
			for (GeneOntologyResult goResult : getProtein().getGoResults()) {
				Pattern p = Pattern.compile(Pattern.quote(token),
						Pattern.CASE_INSENSITIVE);
				if (p.matcher(goResult.getName()).find())
					goScore += goResult.getProbability();
			}
		}
		return goScore;
	}

	public Protein getProtein() {
		return protein;
	}

	public void setProtein(Protein protein) {
		this.protein = protein;
	}

}
