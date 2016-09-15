package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class EvaluationScoreCalculator {

	private Protein protein;
	/**
	 * A ReferenceDescription holds the original Description for Training and
	 * Evaluation-Purposes:
	 */
	private ReferenceDescription referenceDescription;
	/**
	 * A Blast2GoAnnotation is instantiated in order to compare AHRD's
	 * performance with Blast2GOs. Blast2Go assigns multiple descriptions.
	 */
	private Set<Blast2GoAnnot> blast2GoAnnots;
	/**
	 * The unchanged BlastResults are not passed through any filter nor
	 * blacklist, they serve only AHRD-evaluation and training-purposes. Only
	 * the best BlastHit per searched Blast-Database is remembered. Comparison
	 * is based on the Hits' Bit-Scores.
	 */
	private Map<String, BlastResult> unchangedBlastResults = new HashMap<String, BlastResult>();
	private Double evalutionScore;
	private Double evalScoreMinBestCompScore;
	private Double truePositivesRate;
	private Double falsePositivesRate;
	private Double highestPossibleEvaluationScore;
	private Set<GOterm> referenceGoAnnoatations;
	private Double goAnnotationScore;

	public EvaluationScoreCalculator(Protein protein) {
		super();
		setProtein(protein);
	}

	/**
	 * Returns cardinality of intersection between assigned Tokens and reference
	 * Tokens.
	 * 
	 * @param assignedTokens
	 * @param referenceTokens
	 * @return Double - The number of shared Tokens
	 */
	public static Double truePositives(Set<String> assignedTokens,
			Set<String> referenceTokens) {
		double tp = 0.0;
		if (assignedTokens != null && !assignedTokens.isEmpty()) {
			for (String assignedTkn : assignedTokens) {
				if (referenceTokens.contains(assignedTkn))
					tp += 1;
			}
		}
		return tp;
	}

	/**
	 * TPR (True-Positives-Rate) := TP / (TP + FN) TPR is also known as recall.
	 * 
	 * TPR = #shared-tokens / #reference-tokens
	 * 
	 * @param assignedTokens
	 * @param referenceTokens
	 * @return Double - True-Positives-Rate
	 */
	public static Double truePositivesRate(Set<String> assignedTokens,
			Set<String> referenceTokens) {
		return truePositives(assignedTokens, referenceTokens)
				/ referenceTokens.size();
	}

	/**
	 * FPR (False-Positves-Rate) := FP / (FP + TN) FPR is also known as the
	 * fall-out.
	 * 
	 * FPR := #(assigned-tokens not in the reference-tokens) / #(allBlastTokens
	 * without the reference-tokens)
	 * 
	 * @param assignedTokens
	 * @param referenceTokens
	 * @param allBlastTokens
	 * @return Double - False-Positives-Rates
	 */
	public static Double falsePositivesRate(Set<String> assignedTokens,
			Set<String> referenceTokens, Set<String> allBlastTokens) {
		// Count false-positives
		double fp = 0;
		for (String asgnTkn : assignedTokens) {
			if (!referenceTokens.contains(asgnTkn))
				fp += 1;
		}
		// Count all negative tokens:
		double an = allBlastTokens.size();
		for (String blastTkn : allBlastTokens) {
			if (referenceTokens.contains(blastTkn))
				an -= 1;
		}
		// Avoid division by zero:
		return an == 0 ? 0 : fp / an;
	}

	/**
	 * We calculate an annotator's performance as the F-Beta-Score, the harmonic
	 * mean of annotation's precision and recall. The Beta-Parameter is
	 * defaulted to 1.0, but can be set to any real value in the input.yml. Set
	 * to e.g. 2.0 means putting more emphasis on recall than precision. See
	 * http://en.wikipedia.org/wiki/F1_score for details.
	 * 
	 * Note, that the number of total-positives is the number of tokens in the
	 * reference:
	 * 
	 * true-positives (tp) := #shared-tokens
	 * 
	 * precision (pr) := true-positives / #assigned-tokens
	 * 
	 * Recall is identical with the True-Positive-Rate: recall (rc) :=
	 * true-positives / #reference-tokens
	 * 
	 * f-beta-score := (1+beta^2) * (pr * rc) / ((beta^2)*pr + rc)
	 * 
	 * @param assignedTkns
	 *            - Tokens of the Description assigned by AHRD or a
	 *            competitor-method
	 * @param referenceTkns
	 *            - Tokens of the Reference
	 * @return Double - F-Beta-Score
	 */
	public static Double fBetaScore(Set<String> assignedTkns,
			Set<String> referenceTkns) {
		// Validate Reference:
		if (referenceTkns == null || referenceTkns.isEmpty())
			throw new IllegalArgumentException(
					"Cannot calculate F1-Score, got an empty set of Reference-Tokens.");
		// Calculate f-beta-score:
		double fBetaScore = 0.0;
		if (assignedTkns != null && !assignedTkns.isEmpty()) {
			double tp = truePositives(assignedTkns, referenceTkns);
			// Avoid division by zero:
			if (tp > 0.0) {
				double pr = tp / assignedTkns.size();
				double rc = tp / referenceTkns.size();
				// F-Beta-Measure is the harmonic mean of precision and recall
				// weighted by param beta:
				Double bSqr = getSettings().getFMeasureBetaParameter()
						* getSettings().getFMeasureBetaParameter();
				fBetaScore = (1 + bSqr) * (pr * rc) / (bSqr * pr + rc);
			}
		}
		return fBetaScore;
	}

	/**
	 * The unchanged BlastResults are not passed through any filter nor
	 * blacklist, they serve only AHRD-evaluation and training-purposes.
	 * 
	 * @param String
	 *            blastDb
	 * @param BlastResult
	 *            br
	 */
	public void addUnchangedBlastResult(String blastDb, BlastResult br) {
		if (!getUnchangedBlastResults().containsKey(blastDb)
				|| getUnchangedBlastResults().get(blastDb).getBitScore() < br
						.getBitScore()) {
			getUnchangedBlastResults().put(blastDb, br);
		}
	}

	/**
	 * "Germany's Next Top Score" is a show in which AHRD's evaluation-score is
	 * subtracted by the best performing competitor, namely the best unchanged
	 * Blast-Hit.
	 */
	public void assignEvlScrsToCompetitors() {
		if (getReferenceDescription() != null
				&& getReferenceDescription().getDescription() != null) {
			// First Competitor is the Description assigned by AHRD itself:
			if (getProtein().getDescriptionScoreCalculator()
					.getHighestScoringBlastResult() != null) {
				// Generate the set of Evaluation-Tokens from the
				// actually assigned Description, WITHOUT filtering each
				// Token with the BLACKLIST:
				getProtein().getDescriptionScoreCalculator()
						.getHighestScoringBlastResult().tokenizeForEvaluation();
				Set<String> hrdEvlTkns = getProtein()
						.getDescriptionScoreCalculator()
						.getHighestScoringBlastResult().getEvaluationTokens();
				// Calculate the Evaluation-Score as the F-Beta-Score:
				setEvalutionScore(fBetaScore(hrdEvlTkns,
						getReferenceDescription().getTokens()));
				// Enable calculation of the ROC-Curve:
				setTruePositivesRate(truePositivesRate(hrdEvlTkns,
						getReferenceDescription().getTokens()));
				setFalsePositivesRate(falsePositivesRate(hrdEvlTkns,
						getReferenceDescription().getTokens(), getProtein()
								.getTokenScoreCalculator().getTokenScores()
								.keySet()));
			} else {
				// Well, no Description assigned means scores ZERO:
				setEvalutionScore(0.0);
				setTruePositivesRate(0.0);
				setFalsePositivesRate(0.0);
			}
			// Other competitors are the best unchanged BlastHits from all
			// performed Blast-Database-Searches:
			Double bestCompEvlScr = 0.0;
			if (getUnchangedBlastResults().size() > 0) {
				for (String blastDatabase : getUnchangedBlastResults().keySet()) {
					BlastResult cmpt = getUnchangedBlastResults().get(
							blastDatabase);
					if (cmpt != null) {
						// Generate the set of Evaluation-Tokens from the
						// actually assigned Description, WITHOUT filtering each
						// Token with the BLACKLIST:
						cmpt.tokenizeForEvaluation();
						cmpt.setEvaluationScore(fBetaScore(
								cmpt.getEvaluationTokens(),
								getReferenceDescription().getTokens()));
						// Find best performing competitor-method:
						if (cmpt.getEvaluationScore() > bestCompEvlScr)
							bestCompEvlScr = cmpt.getEvaluationScore();
					}
				}
			}
			// Also compare with the Blast2GO-Annotation(s), if present:
			if (getBlast2GoAnnots() != null) {
				for (Blast2GoAnnot b2ga : getBlast2GoAnnots()) {
					b2ga.setEvaluationScore(fBetaScore(
							b2ga.getEvaluationTokens(),
							getReferenceDescription().getTokens()));
					// Find best performing competitor-method:
					if (b2ga.getEvaluationScore() > bestCompEvlScr)
						bestCompEvlScr = b2ga.getEvaluationScore();
				}
			}
			// Compare AHRD's performance:
			setEvalScoreMinBestCompScore(getEvalutionScore() - bestCompEvlScr);
		}
		// Evaluate GO annotations
		setGoAnnotationScore(calcGoAnnotationScore());
	}

	private Double calcGoAnnotationScore() {
		Double f1 = 0.0;
		// TODO Auto-generated method stub
		return f1;
	}

	/**
	 * In order to get more accurate information of how well AHRD performs, we
	 * infer the highest possible score by calculating the evaluation-score for
	 * each BlastResult's Description and remembering the highest achieved
	 * score.
	 */
	public void findHighestPossibleEvaluationScore() {
		setHighestPossibleEvaluationScore(0.0);
		for (List<BlastResult> resultsFromBlastDatabase : getProtein()
				.getBlastResults().values()) {
			for (BlastResult cmpt : resultsFromBlastDatabase) {
				// Generate the set of Evaluation-Tokens from the
				// actually assigned Description, WITHOUT filtering each
				// Token with the BLACKLIST:
				cmpt.tokenizeForEvaluation();
				cmpt.setEvaluationScore(fBetaScore(cmpt.getEvaluationTokens(),
						getReferenceDescription().getTokens()));
				// Find best performing BlastResult-Description:
				if (cmpt.getEvaluationScore() > getHighestPossibleEvaluationScore())
					setHighestPossibleEvaluationScore(cmpt.getEvaluationScore());
			}
		}
	}

	/**
	 * Sorts Blast2GoAnnots by their evaluation-scores ascending. So the
	 * <b>last</b> Blast2GoAnnot in the list will be best performing!
	 * 
	 * @return Sorted List<Blast2GoAnnot>
	 */
	public List<Blast2GoAnnot> sortBlast2GoAnnotsByEvalScore() {
		List<Blast2GoAnnot> b2gaRanked = null;
		if (getBlast2GoAnnots() != null) {
			b2gaRanked = new ArrayList<Blast2GoAnnot>(getBlast2GoAnnots());
			Collections.sort(b2gaRanked);
		}
		return b2gaRanked;
	}

	/**
	 * Blast2Go assigns multiple annotations. In the standard output-file
	 * (.annot) these are written into several lines. Some of those lines show
	 * identical descriptions, hence this method only adds Blast2GoAnnots, if
	 * its description is not already present.
	 * 
	 * @param b2ga
	 */
	public void addBlast2GoAnnot(Blast2GoAnnot b2ga) {
		if (getBlast2GoAnnots() == null)
			setBlast2GoAnnots(new HashSet<Blast2GoAnnot>());
		getBlast2GoAnnots().add(b2ga);
	}

	public ReferenceDescription getReferenceDescription() {
		return referenceDescription;
	}

	public void setReferenceDescription(
			ReferenceDescription referenceDescription) {
		this.referenceDescription = referenceDescription;
	}

	public Map<String, BlastResult> getUnchangedBlastResults() {
		return unchangedBlastResults;
	}

	public void setUnchangedBlastResults(
			Map<String, BlastResult> unchangedBlastResults) {
		this.unchangedBlastResults = unchangedBlastResults;
	}

	public Protein getProtein() {
		return protein;
	}

	public void setProtein(Protein protein) {
		this.protein = protein;
	}

	public Double getEvalutionScore() {
		return evalutionScore;
	}

	public void setEvalutionScore(Double evalutionScore) {
		this.evalutionScore = evalutionScore;
	}

	public Double getEvalScoreMinBestCompScore() {
		return evalScoreMinBestCompScore;
	}

	public void setEvalScoreMinBestCompScore(Double evalScoreMinBestCompScore) {
		this.evalScoreMinBestCompScore = evalScoreMinBestCompScore;
	}

	public Double getTruePositivesRate() {
		return truePositivesRate;
	}

	public void setTruePositivesRate(Double truePositivesRate) {
		this.truePositivesRate = truePositivesRate;
	}

	public Double getFalsePositivesRate() {
		return falsePositivesRate;
	}

	public void setFalsePositivesRate(Double falsePositivesRate) {
		this.falsePositivesRate = falsePositivesRate;
	}

	public Set<Blast2GoAnnot> getBlast2GoAnnots() {
		return blast2GoAnnots;
	}

	public void setBlast2GoAnnots(Set<Blast2GoAnnot> blast2GoAnnots) {
		this.blast2GoAnnots = blast2GoAnnots;
	}

	public Double getHighestPossibleEvaluationScore() {
		return highestPossibleEvaluationScore;
	}

	public void setHighestPossibleEvaluationScore(
			Double highestPossibleEvaluationScore) {
		this.highestPossibleEvaluationScore = highestPossibleEvaluationScore;
	}

	public Set<GOterm> getReferenceGoAnnoatations() {
		return referenceGoAnnoatations;
	}

	public void setReferenceGoAnnoatations(Set<GOterm> referenceGoAnnoatations) {
		this.referenceGoAnnoatations = referenceGoAnnoatations;
	}
	
	public void addReferenceGoAnnotation(GOterm term) {
		if (getReferenceGoAnnoatations() == null)
			setReferenceGoAnnoatations(new HashSet<GOterm>());
		getReferenceGoAnnoatations().add(term);
	}

	public Double getGoAnnotationScore() {
		return goAnnotationScore;
	}

	public void setGoAnnotationScore(Double goAnnotationScore) {
		this.goAnnotationScore = goAnnotationScore;
	}
}
