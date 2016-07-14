package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class DescriptionScoreCalculator {

	private Protein protein;
	private double maxBitScore = 0.0;
	private BlastResult highestScoringBlastResult;
	private Double descriptionHighScore;

	public DescriptionScoreCalculator(Protein protein) {
		super();
		setProtein(protein);
	}

	/**
	 * Assigns each BlastResult's Description-Line its AHRD-Score and then finds
	 * the highest scoring one.
	 * 
	 * @param referenceGoAnnotations
	 *            Map of BlastResults' shortAccesions as keys and their Sets of
	 *            annotated GO Terms as values. If NOT null and any of the query
	 *            proteins' hits of GO Term annotations, AHRD will use the
	 *            highest scoring BlastResult with GO Terms to annotate the
	 *            query.
	 */
	public void findHighestScoringBlastResult(Map<String, Set<String>> referenceGoAnnotations) {
		BlastResult bestScoringBr = null;
		Set<Double> scoreRankingWithGoAnnos = new HashSet<Double>();
		Map<Double, BlastResult> scoreRanking = new HashMap<Double, BlastResult>();
		for (String blastDb : getProtein().getBlastResults().keySet()) {
			for (BlastResult iterBlastResult : getProtein().getBlastResults().get(blastDb)) {
				getProtein().getDescriptionScoreCalculator().calcDescriptionScore(iterBlastResult);
				// Only take Description-Lines into account
				// that have at least a single non-blacklisted Token:
				if (iterBlastResult.getTokens().size() > 0) {
					scoreRanking.put(iterBlastResult.getDescriptionScore(), iterBlastResult);
					if (referenceGoAnnotations != null && !referenceGoAnnotations.isEmpty()
							&& referenceGoAnnotations.containsKey(iterBlastResult.getShortAccession()))
						scoreRankingWithGoAnnos.add(iterBlastResult.getDescriptionScore());
				}
			}
		}
		if (scoreRanking.size() > 0) {
			Set<Double> usedScoreRankin = scoreRankingWithGoAnnos.isEmpty() ? scoreRanking.keySet()
					: scoreRankingWithGoAnnos;
			setDescriptionHighScore(Collections.max(usedScoreRankin));
			bestScoringBr = scoreRanking.get(getDescriptionHighScore());
		}
		setHighestScoringBlastResult(bestScoringBr);
	}

	public void calcDescriptionScore(BlastResult blastResult) {
		blastResult.setDescriptionScore(
				getProtein().getLexicalScoreCalculator().lexicalScore(blastResult) + relativeBlastScore(blastResult));
	}

	public double relativeBlastScore(BlastResult br) {
		return getSettings().getDescriptionScoreBitScoreWeight(br.getBlastDatabaseName()) * br.getBitScore()
				/ getMaxBitScore();
	}

	public void measureMaxBitScore(double bitScore) {
		if (bitScore > getMaxBitScore())
			setMaxBitScore(bitScore);
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

	public void setHighestScoringBlastResult(BlastResult highestScoringBlastResult) {
		this.highestScoringBlastResult = highestScoringBlastResult;
	}

	public Double getDescriptionHighScore() {
		return descriptionHighScore;
	}

	public void setDescriptionHighScore(Double descriptionHighScore) {
		this.descriptionHighScore = descriptionHighScore;
	}
}