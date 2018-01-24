package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class EvaluationScoreCalculator {

	private Protein protein;
	/**
	 * A GroundTruthDescription holds the original Description for Training and
	 * Evaluation-Purposes:
	 */
	private GroundTruthDescription groundTruthDescription;
	/**
	 * The annotations of competitors are instantiated in order to compare AHRD's
	 * performance to them.
	 */
	private Map<String, CompetitorAnnotation> competitorAnnotations;
	
	/**
	 * The unchanged BlastResults are not passed through any filter nor
	 * blacklist, they serve only AHRD-evaluation and training-purposes. Only
	 * the best BlastHit per searched Blast-Database is remembered. Comparison
	 * is based on the Hits' Bit-Scores.
	 */
	private Map<String, BlastResult> bestUnchangedBlastResults = new HashMap<String, BlastResult>();
	private Fscore evalutionScore;
	private Double evalScoreMinBestCompScore;
	private Fscore highestPossibleDescriptionScore;
	private BlastResult blastResultWithHighestPossibleDescriptionScore;
	private Set<GOterm> groundTruthGoAnnoatations = new HashSet<GOterm>();
	private Fscore simpleGoAnnotationScore;
	private Fscore ancestryGoAnnotationScore;
	private Fscore semSimGoAnnotationScore;
	private Fscore highestPossibleSimpleGoAnnotationScore;
	private Fscore highestPossibleAncestryGoAnnotationScore;
	private Fscore highestPossibleSemSimGoAnnotationScore;

	public EvaluationScoreCalculator(Protein protein) {
		super();
		setProtein(protein);
	}

	/**
	 * Returns cardinality of intersection between assigned Tokens and ground truth Tokens.
	 * 
	 * @param assignedTokens
	 * @param groundTruthTokens
	 * @return Double - The number of shared Tokens
	 */
	public static Double truePositives(Set<String> assignedTokens, Set<String> groundTruthTokens) {
		double tp = 0.0;
		if (assignedTokens != null && !assignedTokens.isEmpty()) {
			for (String assignedTkn : assignedTokens) {
				if (groundTruthTokens.contains(assignedTkn))
					tp += 1;
			}
		}
		return tp;
	}

	/**
	 * FPR (False-Positves-Rate) := FP / (FP + TN) FPR is also known as the
	 * fall-out.
	 * 
	 * FPR := #(assigned-tokens not in the ground-truth-tokens) / #(allBlastTokens
	 * without the ground-truth-tokens)
	 * 
	 * @param assignedTokens
	 * @param groundTruthTokens
	 * @param allBlastTokens
	 * @return Double - False-Positives-Rates
	 */
	public static Double falsePositivesRate(Set<String> assignedTokens, Set<String> groundTruthTokens,
			Set<String> allBlastTokens) {
		// Count false-positives
		double fp = 0;
		for (String asgnTkn : assignedTokens) {
			if (!groundTruthTokens.contains(asgnTkn))
				fp += 1;
		}
		// Count all negative tokens:
		double an = allBlastTokens.size();
		for (String blastTkn : allBlastTokens) {
			if (groundTruthTokens.contains(blastTkn))
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
	 * ground truth:
	 * 
	 * true-positives (tp) := #shared-tokens
	 * 
	 * precision (pr) := true-positives / #assigned-tokens
	 * 
	 * Recall is identical with the True-Positive-Rate: recall (rc) :=
	 * true-positives / #ground-truth-tokens
	 * 
	 * f-beta-score := (1+beta^2) * (pr * rc) / ((beta^2)*pr + rc)
	 * 
	 * @param assignedTkns
	 *            - Tokens of the Description assigned by AHRD or a
	 *            competitor-method
	 * @param groundTruthTkns
	 *            - Tokens of the Ground Truth
	 * @return Double - F-Beta-Score
	 */
	public static Fscore fBetaScore(Set<String> assignedTkns, Set<String> groundTruthTkns) {
		// Validate Ground Truth:
		if (groundTruthTkns == null || groundTruthTkns.isEmpty())
			throw new IllegalArgumentException("Cannot calculate F1-Score, got an empty set of Ground-Truth-Tokens.");
		// Calculate f-beta-score:
		Fscore fBetaScore = new Fscore();
		if (assignedTkns != null && !assignedTkns.isEmpty()) {
			double tp = truePositives(assignedTkns, groundTruthTkns);
			// Avoid division by zero:
			if (tp > 0.0) {
				double pr = tp / assignedTkns.size();
				double rc = tp / groundTruthTkns.size();
				fBetaScore.setPrecision(pr);
				fBetaScore.setRecall(rc);
			}
		}
		return fBetaScore;
	}

	/**
	 * The unchanged BlastResults are not passed through any filter nor
	 * blacklist, they serve only AHRD-evaluation and training-purposes.
	 * The BlastResult is compared to the (until now) best BlastResult 
	 * and takes its place if it trumps its bit score. 
	 * 
	 * @param String
	 *            blastDb
	 * @param BlastResult
	 *            br
	 */
	public void addBestUnchangedBlastResult(String blastDb, BlastResult br) {
		if (!getBestUnchangedBlastResults().containsKey(blastDb)
				|| getBestUnchangedBlastResults().get(blastDb).getBitScore() < br.getBitScore()) {
			getBestUnchangedBlastResults().put(blastDb, br);
		}
	}

	/**
	 * "Germany's Next Top Score" is a show in which AHRD's evaluation-score is
	 * subtracted by the best performing competitor, namely the best unchanged
	 * Blast-Hit.
	 */
	public void assignEvaluationScores() {
		if (getGroundTruthDescription() != null && getGroundTruthDescription().getDescription() != null) {
			// First Competitor is the Description assigned by AHRD itself:
			if (getProtein().getDescriptionScoreCalculator().getHighestScoringBlastResult() != null) {
				// Generate the set of evaluation-tokens from the actually assigned description.
				// If evaluateValidTokens is set to false: WITHOUT filtering each token with the BLACKLIST.
				getProtein().getDescriptionScoreCalculator().getHighestScoringBlastResult().tokenizeForEvaluation();
				Set<String> hrdEvlTkns = getProtein().getDescriptionScoreCalculator().getHighestScoringBlastResult()
						.getEvaluationTokens();
				// Calculate the Evaluation-Score as the F-Beta-Score (including Precision and Recall):
				setEvalutionScore(fBetaScore(hrdEvlTkns, getGroundTruthDescription().getTokens()));
			} else {
				// Well, no Description assigned means scores ZERO:
				setEvalutionScore(new Fscore());
			}
			// Do the competitors
			Double bestCompEvlScr = 0.0;
			if (getSettings().hasCompetitors()) {
				for (String competitor : getSettings().getCompetitorSettings().keySet()) {
					Map<String, CompetitorAnnotation> compAnnots = getCompetitorAnnotations();
					if (compAnnots != null) {
						CompetitorAnnotation annot = compAnnots.get(competitor);
						if (annot != null) {
							// Description
							annot.setEvaluationScore(fBetaScore(annot.getEvaluationTokens(), getGroundTruthDescription().getTokens()));
							// Find best performing competitor-method:
							if (annot.getEvaluationScore().getScore() > bestCompEvlScr)
								bestCompEvlScr = annot.getEvaluationScore().getScore();
							// GOterm annotations scores
							if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
								if (getSettings().doCalculateSimpleGoF1Scores())
									annot.setSimpleGoAnnotationScore(calcSimpleGoAnnotationScore(this.groundTruthGoAnnoatations, annot.getGoAnnotations()));
								if (getSettings().doCalculateAncestryGoF1Scores())
									annot.setAncestryGoAnnotationScore(calcAncestryGoAnnotationScore(this.groundTruthGoAnnoatations, annot.getGoAnnotations()));
								if (getSettings().doCalculateSemSimGoF1Scores())
									annot.setSemSimGoAnnotationScore(calcSemSimGoAnnotationScore(this.groundTruthGoAnnoatations, annot.getGoAnnotations()));
							}
						}
					}
				}
			}
			// Other competitors are the best unchanged BlastHits from all Blast-Database-Searches:
			if (getBestUnchangedBlastResults().size() > 0) {
				for (String blastDatabase : getBestUnchangedBlastResults().keySet()) {
					BlastResult cmpt = getBestUnchangedBlastResults().get(blastDatabase);
					if (cmpt != null) {
						// If evaluateValidTokens is set to true: Filter each token with the BLACKLIST.
						if (getSettings().getEvaluateOnlyValidTokens()) {
							cmpt.setEvaluationTokens(TokenScoreCalculator.tokenize(cmpt.getDescription(), getSettings().getTokenBlacklist(cmpt.getBlastDatabaseName())));
						} else {
							cmpt.setEvaluationTokens(cmpt.getTokens());
						}
						cmpt.setEvaluationScore(fBetaScore(cmpt.getEvaluationTokens(), getGroundTruthDescription().getTokens()));
						// Find best performing competitor-method:
						if (cmpt.getEvaluationScore().getScore() > bestCompEvlScr)
							bestCompEvlScr = cmpt.getEvaluationScore().getScore();
						// If requested: Evaluate the GO annotations of the best blast result
						if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
							if (getSettings().doCalculateSimpleGoF1Scores())
								cmpt.setSimpleGoAnnotationScore(
										calcSimpleGoAnnotationScore(this.groundTruthGoAnnoatations, cmpt.getGoAnnotations()));
							if (getSettings().doCalculateAncestryGoF1Scores())
								cmpt.setAncestryGoAnnotationScore(
										calcAncestryGoAnnotationScore(this.groundTruthGoAnnoatations, cmpt.getGoAnnotations()));
							if (getSettings().doCalculateSemSimGoF1Scores())
								cmpt.setSemSimGoAnnotationScore(
										calcSemSimGoAnnotationScore(this.groundTruthGoAnnoatations, cmpt.getGoAnnotations()));
						}
					}
				}
			}
			// Compare AHRD's performance:
			setEvalScoreMinBestCompScore(getEvalutionScore().getScore() - bestCompEvlScr);
		}
		// Evaluate AHRD's GO annotations
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			// Calculation of an F1-score based on ground truth and prediction GO
			// annotations alone
			if (getSettings().doCalculateSimpleGoF1Scores())
				setSimpleGoAnnotationScore(calcSimpleGoAnnotationScore(this.groundTruthGoAnnoatations, this.protein.getGoResultsTerms()));
			// Calculation of an F1-score based on ground truth and prediction GO
			// annotations extended to their complete ancestry
			if (getSettings().doCalculateAncestryGoF1Scores())
				setAncestryGoAnnotationScore(calcAncestryGoAnnotationScore(this.groundTruthGoAnnoatations, this.protein.getGoResultsTerms()));
			// Calculation of an F1-score based on the semantic similarity
			// (based on term information content) of ground truth and prediction
			// GO annotations.
			if (getSettings().doCalculateSemSimGoF1Scores())
				setSemSimGoAnnotationScore(calcSemSimGoAnnotationScore(this.groundTruthGoAnnoatations, this.protein.getGoResultsTerms()));
		}
	}

	/**
	 * Calculates an f score from recall and precision based on simple
	 * cardinality of the ground truth and prediction GO term sets.
	 * 
	 * @return Double - F score
	 */
	private Fscore calcSimpleGoAnnotationScore(Set<GOterm> groundTruth, Set<GOterm> prediction) {
		Fscore f = new Fscore();
		int truePositive = 0;
		if (groundTruth.size() > 0 && prediction.size() > 0) {
			for (Iterator<GOterm> groundTruthIter = groundTruth.iterator(); groundTruthIter.hasNext();) {
				if (prediction.contains(groundTruthIter.next())) {
					truePositive++;
				}
			}
			f.setRecall((double) truePositive / groundTruth.size());
			f.setPrecision((double) truePositive / prediction.size());
		} else {
			if (groundTruth.size() == 0) {
				f.setRecall(1.0);
			}
			if (prediction.size() == 0) {
				f.setPrecision(1.0);
			}
		}
		return f;
	}

	/**
	 * Calculates an f score from recall and precision based on the cardinality
	 * of the ground truth and prediction GO term sets expanded to their complete
	 * ancestry
	 * 
	 * @return Double - F score
	 */
	private Fscore calcAncestryGoAnnotationScore(Set<GOterm> groundTruth, Set<GOterm> prediction) {
		Fscore f = new Fscore();
		Set<GOterm> groundTruthAncestry = new HashSet<GOterm>();
		for (Iterator<GOterm> groundTruthIter = groundTruth.iterator(); groundTruthIter.hasNext();) {
			groundTruthAncestry.addAll(groundTruthIter.next().getAncestry());
		}
		Set<GOterm> predictionAncestry = new HashSet<GOterm>();
		for (Iterator<GOterm> predictionIter = prediction.iterator(); predictionIter.hasNext();) {
			predictionAncestry.addAll(predictionIter.next().getAncestry());
		}
		int truePositive = 0;
		Double recall = 0.0;
		Double precision = 0.0;
		if (groundTruthAncestry.size() > 0 && predictionAncestry.size() > 0) {
			for (Iterator<GOterm> groundTruthIter = groundTruthAncestry.iterator(); groundTruthIter.hasNext();) {
				if (predictionAncestry.contains(groundTruthIter.next())) {
					truePositive++;
				}
			}
			recall = (double) truePositive / groundTruthAncestry.size();
			precision = (double) truePositive / predictionAncestry.size();
		} else {
			if (groundTruthAncestry.size() == 0) {
				recall = 1.0;
			}
			if (predictionAncestry.size() == 0) {
				precision = 1.0;
			}
		}
		f.setPrecision(precision);
		f.setRecall(recall);
		return f;
	}

	/**
	 * Calculates an 1 score from recall and precision based on the semantic
	 * similarity of the ground truth and prediction GO term sets.
	 * 
	 * @return Double - F score
	 */
	private Fscore calcSemSimGoAnnotationScore(Set<GOterm> groundTruth, Set<GOterm> prediction) {
		Fscore f = new Fscore();
		Double recall = 0.0;
		Double precision = 0.0;
		/**
		 * Find the information content of the ground truth: Usually the
		 * information content of the terms themselves, but if they havn't been
		 * used in swissprot their info content is infinite so the highest non
		 * infinite info content from their ancestry is used as fallback.
		 */
		if (groundTruth.size() > 0 && prediction.size() > 0) {
			// recall
			Double infoContentGroundTruth = 0.0;
			Double infoContentPrediction = 0.0;
			for (Iterator<GOterm> groundTruthIter = groundTruth.iterator(); groundTruthIter.hasNext();) {
				GOterm groundTruthTerm = groundTruthIter.next();
				Double maxInfoContentAncestry = 0.0;
				for (Iterator<GOterm> ancestryIter = groundTruthTerm.getAncestry().iterator(); ancestryIter.hasNext();) {
					Double infoContent = ancestryIter.next().getInformationContent();
					if (!Double.isInfinite(infoContent) && infoContent > maxInfoContentAncestry) {
						maxInfoContentAncestry = infoContent;
					}
				}
				infoContentGroundTruth += maxInfoContentAncestry;
				Double maxCommonInfoContentPrediction = 0.0;
				for (Iterator<GOterm> predictionIter = prediction.iterator(); predictionIter
						.hasNext();) {
					Double infoContent = maxCommonInfoContent(groundTruthTerm, predictionIter.next());
					if (!Double.isInfinite(infoContent) && infoContent > maxCommonInfoContentPrediction) {
						maxCommonInfoContentPrediction = infoContent;
					}
				}
				infoContentPrediction += maxCommonInfoContentPrediction;
			}
			if (infoContentGroundTruth > 0.0) {
				recall = infoContentPrediction / infoContentGroundTruth;
			} else {
				recall = 1.0;
			}
			//precision
			infoContentPrediction = 0.0;
			infoContentGroundTruth = 0.0;
			for (Iterator<GOterm> predictionIter = prediction.iterator(); predictionIter
					.hasNext();) {
				GOterm predictionTerm = predictionIter.next();
				Double maxInfoContentAncestry = 0.0;
				for (Iterator<GOterm> ancestryIter = predictionTerm.getAncestry().iterator(); ancestryIter.hasNext();) {
					Double infoContent = ancestryIter.next().getInformationContent();
					if (!Double.isInfinite(infoContent) && infoContent > maxInfoContentAncestry) {
						maxInfoContentAncestry = infoContent;
					}
				}
				infoContentPrediction += maxInfoContentAncestry;
				Double maxCommonInfoContentGroundTruth = 0.0;
				for (Iterator<GOterm> groundTruthIter = groundTruth.iterator(); groundTruthIter
						.hasNext();) {
					Double infoContent = maxCommonInfoContent(predictionTerm, groundTruthIter.next());
					if (!Double.isInfinite(infoContent) && infoContent > maxCommonInfoContentGroundTruth) {
						maxCommonInfoContentGroundTruth = infoContent;
					}
				}
				infoContentGroundTruth += maxCommonInfoContentGroundTruth;
			}
			if (infoContentPrediction > 0.0) {
				precision = infoContentGroundTruth / infoContentPrediction;
			} else {
				precision = 1.0;
			}
		} else {
			if (groundTruth.size() == 0) {
				recall = 1.0;
			}
			if (prediction.size() == 0) {
				precision = 1.0;
			}
		}
		if (f.getScore() > 1.0) { // Something went very wrong - Should never be happening
			System.out.println("p: " + precision + "\tr: " + recall + "\tf1: " + f);
		}
		f.setPrecision(precision);
		f.setRecall(recall);
		return f;
	}

	private Double maxCommonInfoContent(GOterm firstTerm, GOterm secondTerm) {
		Set<GOterm> commonAncestry = new HashSet<GOterm>(firstTerm.getAncestry());
		commonAncestry.retainAll(secondTerm.getAncestry());
		Double maxCommonInfoContent = 0.0;
		for (Iterator<GOterm> ancestryIter = commonAncestry.iterator(); ancestryIter.hasNext();) {
			Double infoContent = ancestryIter.next().getInformationContent();
			if (!Double.isInfinite(infoContent) && infoContent > maxCommonInfoContent) {
				maxCommonInfoContent = infoContent;
			}
		}
		return maxCommonInfoContent;
	}

	/**
	 * In order to get more accurate information of how well AHRD performs, we
	 * infer the highest possible score by calculating the evaluation-score for
	 * each BlastResult's Description and remembering the highest achieved
	 * score.
	 */
	public void findBlastResultWithHighestPossibleDescriptionScore() {
		setHighestPossibleDescriptionScore(new Fscore());
		for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
			for (BlastResult cmpt : resultsFromBlastDatabase) {
				// Generate the set of evaluation-tokens for each description, 
				// if evaluateValidTokens is set to false, WITHOUT filtering each token with the BLACKLIST.
				cmpt.tokenizeForEvaluation();
				cmpt.setEvaluationScore(fBetaScore(cmpt.getEvaluationTokens(), getGroundTruthDescription().getTokens()));
				// Find best performing BlastResult-Description:
				if (cmpt.getEvaluationScore().getScore() > getHighestPossibleDescriptionScore().getScore()) {
					setHighestPossibleDescriptionScore(cmpt.getEvaluationScore());
					setBlastResultWithHighestPossibleDescriptionScore(cmpt);
				}
			}
		}
	}
	
	/**
	 * In order to get more accurate information of how well AHRD performs, we
	 * infer the highest possible score by calculating all requested 
	 * GO-annotation-scores (simple, ancestry, semsim) for each 
	 * BlastResult's GO annotations and remembering the highest achieved score.
	 */
	public void findHighestPossibleGoScore() {
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			if (getSettings().doCalculateSimpleGoF1Scores()) {
				this.setHighestPossibleSimpleGoAnnotationScore(calcSimpleGoAnnotationScore(this.groundTruthGoAnnoatations, new HashSet<GOterm>())); // In case ground truth go annotation is empty
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcSimpleGoAnnotationScore(this.groundTruthGoAnnoatations, br.getGoAnnotations());
						if (score.getScore() > this.getHighestPossibleSimpleGoAnnotationScore().getScore())
								this.setHighestPossibleSimpleGoAnnotationScore(score);
						if (getHighestPossibleSimpleGoAnnotationScore().getScore().equals(1.0))
							break;
					}
					if (getHighestPossibleSimpleGoAnnotationScore().getScore().equals(1.0))
						break;
				}
			}
			if (getSettings().doCalculateAncestryGoF1Scores()) {
				this.setHighestPossibleAncestryGoAnnotationScore(this.calcAncestryGoAnnotationScore(this.groundTruthGoAnnoatations, new HashSet<GOterm>())); // In case ground truth go annotation is empty
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcAncestryGoAnnotationScore(this.groundTruthGoAnnoatations, br.getGoAnnotations());
						if (score.getScore() > this.getHighestPossibleAncestryGoAnnotationScore().getScore())
								this.setHighestPossibleAncestryGoAnnotationScore(score);
						if (getHighestPossibleAncestryGoAnnotationScore().getScore().equals(1.0))
							break;
					}
					if (getHighestPossibleAncestryGoAnnotationScore().getScore().equals(1.0))
						break;
				}
			}
			if (getSettings().doCalculateSemSimGoF1Scores()) {
				this.setHighestPossibleSemSimGoAnnotationScore(calcSemSimGoAnnotationScore(this.groundTruthGoAnnoatations, new HashSet<GOterm>())); // In case ground truth go annotation is empty
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcSemSimGoAnnotationScore(this.groundTruthGoAnnoatations, br.getGoAnnotations());
						if (score.getScore() > this.getHighestPossibleSemSimGoAnnotationScore().getScore())
								this.setHighestPossibleSemSimGoAnnotationScore(score);
						if (getHighestPossibleSemSimGoAnnotationScore().getScore().equals(1.0))
							break;
					}
					if (getHighestPossibleSemSimGoAnnotationScore().getScore().equals(1.0))
						break;
				}
			}
		}
	}

	public GroundTruthDescription getGroundTruthDescription() {
		return groundTruthDescription;
	}

	public void setGroundTruthDescription(GroundTruthDescription groundTruthDescription) {
		this.groundTruthDescription = groundTruthDescription;
	}

	public Map<String, BlastResult> getBestUnchangedBlastResults() {
		return bestUnchangedBlastResults;
	}

	public void setBestUnchangedBlastResults(Map<String, BlastResult> unchangedBlastResults) {
		this.bestUnchangedBlastResults = unchangedBlastResults;
	}

	public Protein getProtein() {
		return protein;
	}

	public void setProtein(Protein protein) {
		this.protein = protein;
	}

	public Fscore getEvalutionScore() {
		return evalutionScore;
	}

	public void setEvalutionScore(Fscore evalutionScore) {
		this.evalutionScore = evalutionScore;
	}

	public Double getEvalScoreMinBestCompScore() {
		return evalScoreMinBestCompScore;
	}

	public void setEvalScoreMinBestCompScore(Double evalScoreMinBestCompScore) {
		this.evalScoreMinBestCompScore = evalScoreMinBestCompScore;
	}

	public Fscore getHighestPossibleDescriptionScore() {
		return highestPossibleDescriptionScore;
	}

	public void setHighestPossibleDescriptionScore(Fscore highestPossibleEvaluationScore) {
		this.highestPossibleDescriptionScore = highestPossibleEvaluationScore;
	}

	public Set<GOterm> getGroundTruthGoAnnoatations() {
		return groundTruthGoAnnoatations;
	}

	public void setGroundTruthGoAnnoatations(Set<GOterm> groundTruthGoAnnoatations) {
		this.groundTruthGoAnnoatations = groundTruthGoAnnoatations;
	}

	public Fscore getSimpleGoAnnotationScore() {
		return simpleGoAnnotationScore;
	}

	public void setSimpleGoAnnotationScore(Fscore simpleGoAnnotationScore) {
		this.simpleGoAnnotationScore = simpleGoAnnotationScore;
	}

	public Fscore getAncestryGoAnnotationScore() {
		return ancestryGoAnnotationScore;
	}

	public void setAncestryGoAnnotationScore(Fscore ancestryGoAnnotationScore) {
		this.ancestryGoAnnotationScore = ancestryGoAnnotationScore;
	}

	public Fscore getSemSimGoAnnotationScore() {
		return semSimGoAnnotationScore;
	}

	public void setSemSimGoAnnotationScore(Fscore semSimGoAnnotationScore) {
		this.semSimGoAnnotationScore = semSimGoAnnotationScore;
	}

	public Map<String, CompetitorAnnotation> getCompetitorAnnotations() {
		return competitorAnnotations;
	}

	public void setCompetitorAnnotations(Map<String, CompetitorAnnotation> competitorAnnotations) {
		this.competitorAnnotations = competitorAnnotations;
	}
	
	public void addCompetitorAnnotation(String competitor, CompetitorAnnotation annotation) {
		if (this.competitorAnnotations == null) {
			setCompetitorAnnotations(new HashMap<String, CompetitorAnnotation>());
		}
		this.competitorAnnotations.put(competitor, annotation);
	}

	public Fscore getHighestPossibleSimpleGoAnnotationScore() {
		return highestPossibleSimpleGoAnnotationScore;
	}

	public void setHighestPossibleSimpleGoAnnotationScore(Fscore highestPossiblesimpleGoAnnotationScore) {
		this.highestPossibleSimpleGoAnnotationScore = highestPossiblesimpleGoAnnotationScore;
	}

	public Fscore getHighestPossibleAncestryGoAnnotationScore() {
		return highestPossibleAncestryGoAnnotationScore;
	}

	public void setHighestPossibleAncestryGoAnnotationScore(Fscore highestPossibleancestryGoAnnotationScore) {
		this.highestPossibleAncestryGoAnnotationScore = highestPossibleancestryGoAnnotationScore;
	}

	public Fscore getHighestPossibleSemSimGoAnnotationScore() {
		return highestPossibleSemSimGoAnnotationScore;
	}

	public void setHighestPossibleSemSimGoAnnotationScore(Fscore highestPossiblesemSimGoAnnotationScore) {
		this.highestPossibleSemSimGoAnnotationScore = highestPossiblesemSimGoAnnotationScore;
	}

	public BlastResult getBlastResultWithHighestPossibleDescriptionScore() {
		return blastResultWithHighestPossibleDescriptionScore;
	}

	public void setBlastResultWithHighestPossibleDescriptionScore(
			BlastResult blastResultWithHighestPossibleDescriptionScore) {
		this.blastResultWithHighestPossibleDescriptionScore = blastResultWithHighestPossibleDescriptionScore;
	}

}
