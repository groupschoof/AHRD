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
			for (BlastResult iterBlastResult : getProtein().getBlastResults().get(blastDb)) {
				getProtein().getDescriptionScoreCalculator().calcDescriptionScore(iterBlastResult);
				// Only take Description-Lines into account
				// that have at least a single non-blacklisted Token:
				if (iterBlastResult.getTokens().size() > 0) {
//					System.out.println(iterBlastResult.getAccession() + "\t"
//							+ TokenScoreCalculator.overlapScore(iterBlastResult.getQueryStart(),
//									iterBlastResult.getQueryEnd(), getProtein().getSequenceLength(),
//									iterBlastResult.getSubjectStart(), iterBlastResult.getSubjectEnd(),
//									iterBlastResult.getSubjectLength())
//							+ "\t" + iterBlastResult.getDescriptionScore());
					scoreRanking.put(iterBlastResult.getDescriptionScore(), iterBlastResult);
				}
			}
		}
		if (scoreRanking.size() > 0) {
			setDescriptionHighScore(Collections.max(scoreRanking.keySet()));
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