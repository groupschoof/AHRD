package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class DescriptionScoreCalculator {

	private Protein protein;
	private double maxBitScore = 0.0;
	private BlastResult highestScoringBlastResult;
	private Double descriptionHighScore;
	private Map<String, Integer> descLinePatternFrequencies = new HashMap<String, Integer>();
	private int maxDescriptionLineFrequency = 0;

	public DescriptionScoreCalculator(Protein protein) {
		super();
		setProtein(protein);
	}

	/**
	 * Assigns each BlastResult's Description-Line its AHRD-Score and then finds
	 * the highest scoring one.
	 */
	public void findHighestScoringBlastResult() {
		BlastResult bestScoringBr = null;
		Map<Double, BlastResult> scoreRanking = new HashMap<Double, BlastResult>();
		for (String blastDb : getProtein().getBlastResults().keySet()) {
			for (BlastResult iterBlastResult : getProtein().getBlastResults()
					.get(blastDb)) {
				getProtein().getDescriptionScoreCalculator()
						.calcDescriptionScore(iterBlastResult);
				// Only take Description-Lines into account
				// that have at least a single non-blacklisted Token:
				if (iterBlastResult.getTokens().size() > 0)
					scoreRanking.put(iterBlastResult.getDescriptionScore(),
							iterBlastResult);
			}
		}
		if (scoreRanking.size() > 0) {
			setDescriptionHighScore(Collections.max(scoreRanking.keySet()));
			bestScoringBr = scoreRanking.get(getDescriptionHighScore());
		}
		setHighestScoringBlastResult(bestScoringBr);
	}

	public void calcDescriptionScore(BlastResult blastResult) {
		blastResult.setDescriptionScore(getProtein()
				.getLexicalScoreCalculator().lexicalScore(blastResult)
				+ relativeBlastScore(blastResult) + patternFactor(blastResult));
	}

	/**
	 * Calculates the PatternFactor, which is the frequency of the description
	 * divided by the number of occurrences of the most frequent description.
	 * 
	 * @param br
	 * @return
	 */
	public double patternFactor(BlastResult br) {
		return getSettings().getDescriptionScorePatternFactorWeight()
				* (getDescLinePatternFrequencies().get(br.patternize()) / getMaxDescriptionLineFrequency());
	}

	public double relativeBlastScore(BlastResult br) {
		return getSettings().getDescriptionScoreBitScoreWeight(
				br.getBlastDatabaseName())
				* br.getBitScore() / getMaxBitScore();
	}

	/**
	 * Computes argument BlastResult's 'domain similarity score' times the
	 * configurable weight. Always returns a valid real number, even if the
	 * BlastResult does not have a domain similarity score assigned. The latter
	 * is the case if there are no Interpro annotations available for either the
	 * query protein or the BlastResult itself or AHRD was launched without any
	 * Interpro annotations at all.
	 * 
	 * @param BlastResult
	 *            br
	 * @return Double
	 */
	public Double domainSimilarityScore(BlastResult br) {
		// (1) Instantiate a dss with value zero.
		// (2) if and only if the BlastResult has a non null domain similarity score
		// multiply it with the configurable weight
		// 'descriptionScoreDomainSimilarityWeight'
		// return result
		return null;
	}

	public void measureMaxBitScore(double bitScore) {
		if (bitScore > getMaxBitScore())
			setMaxBitScore(bitScore);
	}

	public void measureDescriptionLineFrequency(String descLinePattern) {
		Integer currentFrequency = getDescLinePatternFrequencies().containsKey(
				descLinePattern) ? getDescLinePatternFrequencies().get(
				descLinePattern) + 1 : 1;
		getDescLinePatternFrequencies().put(descLinePattern, currentFrequency);
		// Measure maxDescriptionLineFrequency:
		if (currentFrequency > getMaxDescriptionLineFrequency())
			setMaxDescriptionLineFrequency(currentFrequency);
	}

	/**
	 * Get protein.
	 * 
	 * @return protein as Protein.
	 */
	public Protein getProtein() {
		return protein;
	}

	/**
	 * Set protein.
	 * 
	 * @param protein
	 *            the value to set.
	 */
	public void setProtein(Protein protein) {
		this.protein = protein;
	}

	/**
	 * Get maxBitScore.
	 * 
	 * @return maxBitScore as double.
	 */
	public double getMaxBitScore() {
		return maxBitScore;
	}

	/**
	 * Set maxBitScore.
	 * 
	 * @param maxBitScore
	 *            the value to set.
	 */
	public void setMaxBitScore(double maxBitScore) {
		this.maxBitScore = maxBitScore;
	}

	public BlastResult getHighestScoringBlastResult() {
		return highestScoringBlastResult;
	}

	public void setHighestScoringBlastResult(
			BlastResult highestScoringBlastResult) {
		this.highestScoringBlastResult = highestScoringBlastResult;
	}

	public Double getDescriptionHighScore() {
		return descriptionHighScore;
	}

	public void setDescriptionHighScore(Double descriptionHighScore) {
		this.descriptionHighScore = descriptionHighScore;
	}

	public Map<String, Integer> getDescLinePatternFrequencies() {
		return descLinePatternFrequencies;
	}

	public void setDescLinePatternFrequencies(
			Map<String, Integer> descriptionLineFrequencies) {
		this.descLinePatternFrequencies = descriptionLineFrequencies;
	}

	/**
	 * Get maxDescriptionLineFrequency.
	 * 
	 * @return maxDescriptionLineFrequency as int.
	 */
	public int getMaxDescriptionLineFrequency() {
		return maxDescriptionLineFrequency;
	}

	/**
	 * Set maxDescriptionLineFrequency.
	 * 
	 * @param maxDescriptionLineFrequency
	 *            the value to set.
	 */
	public void setMaxDescriptionLineFrequency(int maxDescriptionLineFrequency) {
		this.maxDescriptionLineFrequency = maxDescriptionLineFrequency;
	}

}
