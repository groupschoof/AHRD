package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import static ahrd.model.AhrdDb.getReferenceProteinDAO;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DescriptionScoreCalculator {

	/**
	 * Global implementation of the Description Blacklist.
	 * 
	 * @param description
	 * @param blacklist
	 * @return TRUE if and only if none of the regular expressions in blacklist
	 *         matches the argument description. FALSE otherwise.
	 */
	public static boolean passesBlacklist(String description, List<String> blacklist) {
		boolean passesBlacklist = (description != null && !description.equals(""));
		for (Iterator<String> i = blacklist.iterator(); (i.hasNext() && passesBlacklist);) {
			Pattern p = Pattern.compile(i.next());
			Matcher m = p.matcher(description);
			passesBlacklist = !m.find();
		}
		return passesBlacklist;
	}

	/**
	 * Global implementation of the filter Description function.
	 * 
	 * @param description
	 * @param filter
	 * @return A modified version of argument description in which all matches
	 *         to any of the regular expressions in argument filter are deleted.
	 *         Finally the filtered description is trimmed and multiple
	 *         white-spaces are condensed into a single white-spaces.
	 */
	public static String filter(String description, List<String> filter) {
		String filteredDescLine = description;
		for (Iterator<String> i = filter.iterator(); i.hasNext();) {
			Pattern p = Pattern.compile(i.next());
			// Replace with whitespace, so word-boundaries are kept up
			filteredDescLine = p.matcher(filteredDescLine).replaceAll(" ");
		}
		// Condense multiple whitespaces into one and trim the description-line:
		filteredDescLine = filteredDescLine.replaceAll("\\s{2,}", " ").trim();
		return filteredDescLine;
	}

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
		ReferenceProtein rp;
		Set<Double> scoreRankingWithGoAnnos = new HashSet<Double>();
		Map<Double, BlastResult> scoreRanking = new HashMap<Double, BlastResult>();
		for (String blastDb : getProtein().getBlastResults().keySet()) {
			for (BlastResult iterBlastResult : getProtein().getBlastResults().get(blastDb)) {
				getProtein().getDescriptionScoreCalculator().calcDescriptionScore(iterBlastResult);
				// Only take Description-Lines into account
				// that have at least a single non-blacklisted Token:
				if (iterBlastResult.getTokens().size() > 0) {
					scoreRanking.put(iterBlastResult.getDescriptionScore(), iterBlastResult);
					rp = getReferenceProteinDAO().byAccession.get(iterBlastResult.getAccession());
					if (rp != null && !rp.getGoTerms().isEmpty() && getSettings().getPreferReferenceWithGoAnnos())
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
		// If AHRD is able to annotate Gene Ontology (GO) Terms, do so:
		if (bestScoringBr != null) {
			rp = getReferenceProteinDAO().byAccession.get(bestScoringBr.getAccession());
			if (rp != null && !rp.getGoTerms().isEmpty())
				getProtein().setGoResults(rp.getGoTerms());
		}
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
