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
	 * A ReferenceDescription holds the original Description for Training and
	 * Evaluation-Purposes:
	 */
	private ReferenceDescription referenceDescription;
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
	private Set<GOterm> referenceGoAnnoatations = new HashSet<GOterm>();
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
	 * Returns cardinality of intersection between assigned Tokens and reference
	 * Tokens.
	 * 
	 * @param assignedTokens
	 * @param referenceTokens
	 * @return Double - The number of shared Tokens
	 */
	public static Double truePositives(Set<String> assignedTokens, Set<String> referenceTokens) {
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
	public static Double falsePositivesRate(Set<String> assignedTokens, Set<String> referenceTokens,
			Set<String> allBlastTokens) {
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
	public static Fscore fBetaScore(Set<String> assignedTkns, Set<String> referenceTkns) {
		// Validate Reference:
		if (referenceTkns == null || referenceTkns.isEmpty())
			throw new IllegalArgumentException("Cannot calculate F1-Score, got an empty set of Reference-Tokens.");
		// Calculate f-beta-score:
		Fscore fBetaScore = new Fscore();
		if (assignedTkns != null && !assignedTkns.isEmpty()) {
			double tp = truePositives(assignedTkns, referenceTkns);
			// Avoid division by zero:
			if (tp > 0.0) {
				double pr = tp / assignedTkns.size();
				double rc = tp / referenceTkns.size();
				// F-Beta-Measure is the harmonic mean of precision and recall
				// weighted by param beta:
				Double bSqr = getSettings().getFMeasureBetaParameter() * getSettings().getFMeasureBetaParameter();
				fBetaScore.setScore((1 + bSqr) * (pr * rc) / (bSqr * pr + rc));
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
	public void assignEvlScrsToCompetitors() {
		if (getReferenceDescription() != null && getReferenceDescription().getDescription() != null) {
			// First Competitor is the Description assigned by AHRD itself:
			if (getProtein().getDescriptionScoreCalculator().getHighestScoringBlastResult() != null) {
				// Generate the set of Evaluation-Tokens from the
				// actually assigned Description, WITHOUT filtering each
				// Token with the BLACKLIST:
				getProtein().getDescriptionScoreCalculator().getHighestScoringBlastResult().tokenizeForEvaluation();
				Set<String> hrdEvlTkns = getProtein().getDescriptionScoreCalculator().getHighestScoringBlastResult()
						.getEvaluationTokens();
				// Calculate the Evaluation-Score as the F-Beta-Score (including Precision and Recall):
				setEvalutionScore(fBetaScore(hrdEvlTkns, getReferenceDescription().getTokens()));
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
							annot.setEvaluationScore(fBetaScore(annot.getEvaluationTokens(), getReferenceDescription().getTokens()));
							// Find best performing competitor-method:
							if (annot.getEvaluationScore().getScore() > bestCompEvlScr)
								bestCompEvlScr = annot.getEvaluationScore().getScore();
							// GOterm annotations scores
							if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
								if (getSettings().doCalculateSimpleGoF1Scores())
									annot.setSimpleGoAnnotationScore(calcSimpleGoAnnotationScore(this.referenceGoAnnoatations, annot.getGoAnnotations()));
								if (getSettings().doCalculateAncestryGoF1Scores())
									annot.setAncestryGoAnnotationScore(calcAncestryGoAnnotationScore(this.referenceGoAnnoatations, annot.getGoAnnotations()));
								if (getSettings().doCalculateSemSimGoF1Scores())
									annot.setSemSimGoAnnotationScore(calcSemSimGoAnnotationScore(this.referenceGoAnnoatations, annot.getGoAnnotations()));
							}
						}
					}
				}
			}
			if (this.protein.getAccession().equals("Q9UTR4")) {
				System.out.println("Whe're here!");
			}
			// Other competitors are the best unchanged BlastHits from all
			// performed Blast-Database-Searches:
			if (getBestUnchangedBlastResults().size() > 0) {
				for (String blastDatabase : getBestUnchangedBlastResults().keySet()) {
					BlastResult cmpt = getBestUnchangedBlastResults().get(blastDatabase);
					if (cmpt != null) {
						// Generate the set of Evaluation-Tokens from the
						// actually assigned Description, WITHOUT filtering each
						// Token with the BLACKLIST:
						cmpt.tokenizeForEvaluation();
						cmpt.setEvaluationScore(
								fBetaScore(cmpt.getEvaluationTokens(), getReferenceDescription().getTokens()));
						// Find best performing competitor-method:
						if (cmpt.getEvaluationScore().getScore() > bestCompEvlScr)
							bestCompEvlScr = cmpt.getEvaluationScore().getScore();
						// If requested: Evaluate the GO annotations of the best blast result
						if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
							if (getSettings().doCalculateSimpleGoF1Scores())
								cmpt.setSimpleGoAnnotationScore(
										calcSimpleGoAnnotationScore(this.referenceGoAnnoatations, cmpt.getGoAnnotations()));
							if (getSettings().doCalculateAncestryGoF1Scores())
								cmpt.setAncestryGoAnnotationScore(
										calcAncestryGoAnnotationScore(this.referenceGoAnnoatations, cmpt.getGoAnnotations()));
							if (getSettings().doCalculateSemSimGoF1Scores())
								cmpt.setSemSimGoAnnotationScore(
										calcSemSimGoAnnotationScore(this.referenceGoAnnoatations, cmpt.getGoAnnotations()));
						}
					}
				}
			}
			// Compare AHRD's performance:
			setEvalScoreMinBestCompScore(getEvalutionScore().getScore() - bestCompEvlScr);
		}
		// Evaluate GO annotations
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
			/**
			 * AHRD
			 */
			// Calculation of an F1-score based on reference and prediction GO
			// annotations alone
			if (getSettings().doCalculateSimpleGoF1Scores())
				setSimpleGoAnnotationScore(calcSimpleGoAnnotationScore(this.referenceGoAnnoatations, this.protein.getGoResultsTerms()));
			// Calculation of an F1-score based on reference and prediction GO
			// annotations extended to their complete ancestry
			if (getSettings().doCalculateAncestryGoF1Scores())
				setAncestryGoAnnotationScore(calcAncestryGoAnnotationScore(this.referenceGoAnnoatations, this.protein.getGoResultsTerms()));
			// Calculation of an F1-score based on the semantic similarity
			// (based on term information content) of reference and prediction
			// GO annotations.
			if (getSettings().doCalculateSemSimGoF1Scores())
				setSemSimGoAnnotationScore(calcSemSimGoAnnotationScore(this.referenceGoAnnoatations, this.protein.getGoResultsTerms()));
		}
	}

	/**
	 * Calculates an f score from recall and precision based on simple
	 * cardinality of the reference and prediction GO term sets.
	 * 
	 * @return Double - F score
	 */
	private Fscore calcSimpleGoAnnotationScore(Set<GOterm> reference, Set<GOterm> prediction) {
		Fscore f = new Fscore();
		int truePositive = 0;
		if (reference.size() > 0 && prediction.size() > 0) {
			for (Iterator<GOterm> referenceIter = reference.iterator(); referenceIter.hasNext();) {
				if (prediction.contains(referenceIter.next())) {
					truePositive++;
				}
			}
			f.setRecall((double) truePositive / reference.size());
			f.setPrecision((double) truePositive / prediction.size());
		} else {
			if (reference.size() == 0) {
				f.setRecall(1.0);
			}
			if (prediction.size() == 0) {
				f.setPrecision(1.0);
			}
		}
		if (f.getPrecision() > 0.0 && f.getRecall() > 0.0) {
			Double bSqr = getSettings().getFMeasureBetaParameter() * getSettings().getFMeasureBetaParameter();
			f.setScore((1 + bSqr) * f.getPrecision() * f.getRecall() / (bSqr * f.getPrecision() + f.getRecall()));
		}
		return f;
	}

	/**
	 * Calculates an f score from recall and precision based on the cardinality
	 * of the reference and prediction GO term sets expanded to their complete
	 * ancestry
	 * 
	 * @return Double - F score
	 */
	private Fscore calcAncestryGoAnnotationScore(Set<GOterm> reference, Set<GOterm> prediction) {
		Fscore f = new Fscore();
		Set<GOterm> referenceAncestry = new HashSet<GOterm>();
		for (Iterator<GOterm> referenceIter = reference.iterator(); referenceIter.hasNext();) {
			referenceAncestry.addAll(referenceIter.next().getAncestry());
		}
		Set<GOterm> predictionAncestry = new HashSet<GOterm>();
		for (Iterator<GOterm> predictionIter = prediction.iterator(); predictionIter.hasNext();) {
			predictionAncestry.addAll(predictionIter.next().getAncestry());
		}
		int truePositive = 0;
		Double recall = 0.0;
		Double precision = 0.0;
		if (referenceAncestry.size() > 0 && predictionAncestry.size() > 0) {
			for (Iterator<GOterm> referenceIter = referenceAncestry.iterator(); referenceIter.hasNext();) {
				if (predictionAncestry.contains(referenceIter.next())) {
					truePositive++;
				}
			}
			recall = (double) truePositive / referenceAncestry.size();
			precision = (double) truePositive / predictionAncestry.size();
		} else {
			if (referenceAncestry.size() == 0) {
				recall = 1.0;
			}
			if (predictionAncestry.size() == 0) {
				precision = 1.0;
			}
		}
		if (precision > 0.0 && recall > 0.0) {
			Double bSqr = getSettings().getFMeasureBetaParameter() * getSettings().getFMeasureBetaParameter();
			f.setScore((1 + bSqr) * precision * recall / (bSqr * precision + recall));
		}
		f.setPrecision(precision);
		f.setRecall(recall);
		return f;
	}

	/**
	 * Calculates an 1 score from recall and precision based on the semantic
	 * similarity of the reference and prediction GO term sets.
	 * 
	 * @return Double - F score
	 */
	private Fscore calcSemSimGoAnnotationScore(Set<GOterm> reference, Set<GOterm> prediction) {
		Fscore f = new Fscore();
		Double recall = 0.0;
		Double precision = 0.0;
		// Recall
		/**
		 * Find the information content of the reference: Usually the
		 * information content of the terms themselves, but if they havn't been
		 * used in swissprot their info content is infinite so the highest non
		 * infinite info content from their ancestry is used as fallback.
		 */
		if (reference.size() > 0 && prediction.size() > 0) {
			Double infoContentReference = 0.0;
			Double infoContentPrediction = 0.0;
			for (Iterator<GOterm> referenceIter = reference.iterator(); referenceIter.hasNext();) {
				GOterm referenceTerm = referenceIter.next();
				Double maxInfoContentAncestry = 0.0;
				for (Iterator<GOterm> ancestryIter = referenceTerm.getAncestry().iterator(); ancestryIter.hasNext();) {
					Double infoContent = ancestryIter.next().getInformationContent();
					if (!Double.isInfinite(infoContent) && infoContent > maxInfoContentAncestry) {
						maxInfoContentAncestry = infoContent;
					}
				}
				infoContentReference += maxInfoContentAncestry;
				Double maxCommonInfoContentPrediction = 0.0;
				for (Iterator<GOterm> predictionIter = prediction.iterator(); predictionIter
						.hasNext();) {
					Double infoContent = maxCommonInfoContent(referenceTerm, predictionIter.next());
					if (!Double.isInfinite(infoContent) && infoContent > maxCommonInfoContentPrediction) {
						maxCommonInfoContentPrediction = infoContent;
					}
				}
				infoContentPrediction += maxCommonInfoContentPrediction;
			}
			if (infoContentReference > 0.0) {
				recall = infoContentPrediction / infoContentReference;
			} else {
				recall = 1.0;
			}

			infoContentPrediction = 0.0;
			infoContentReference = 0.0;
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
				Double maxCommonInfoContentReference = 0.0;
				for (Iterator<GOterm> referenceIter = reference.iterator(); referenceIter
						.hasNext();) {
					Double infoContent = maxCommonInfoContent(predictionTerm, referenceIter.next());
					if (!Double.isInfinite(infoContent) && infoContent > maxCommonInfoContentReference) {
						maxCommonInfoContentReference = infoContent;
					}
				}
				infoContentReference += maxCommonInfoContentReference;

			}
			if (infoContentPrediction > 0.0) {
				precision = infoContentReference / infoContentPrediction;
			} else {
				precision = 1.0;
			}
		} else {
			if (reference.size() == 0) {
				recall = 1.0;
			}
			if (prediction.size() == 0) {
				precision = 1.0;
			}
		}
		if (precision > 0.0 && recall > 0.0) {
			Double bSqr = getSettings().getFMeasureBetaParameter() * getSettings().getFMeasureBetaParameter();
			f.setScore((1 + bSqr) * precision * recall / (bSqr * precision + recall));
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
				// Generate the set of Evaluation-Tokens from the
				// actually assigned Description, WITHOUT filtering each
				// Token with the BLACKLIST:
				cmpt.tokenizeForEvaluation();
				cmpt.setEvaluationScore(fBetaScore(cmpt.getEvaluationTokens(), getReferenceDescription().getTokens()));
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
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
			if (getSettings().doCalculateSimpleGoF1Scores()) {
				this.setHighestPossibleSimpleGoAnnotationScore(calcSimpleGoAnnotationScore(this.referenceGoAnnoatations, new HashSet<GOterm>())); // In case reference go annotation is empty
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcSimpleGoAnnotationScore(this.referenceGoAnnoatations, br.getGoAnnotations());
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
				this.setHighestPossibleAncestryGoAnnotationScore(this.calcAncestryGoAnnotationScore(this.referenceGoAnnoatations, new HashSet<GOterm>())); // In case reference go annotation is empty
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcAncestryGoAnnotationScore(this.referenceGoAnnoatations, br.getGoAnnotations());
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
				this.setHighestPossibleSemSimGoAnnotationScore(calcSemSimGoAnnotationScore(this.referenceGoAnnoatations, new HashSet<GOterm>())); // In case reference go annotation is empty
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcSemSimGoAnnotationScore(this.referenceGoAnnoatations, br.getGoAnnotations());
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

	public ReferenceDescription getReferenceDescription() {
		return referenceDescription;
	}

	public void setReferenceDescription(ReferenceDescription referenceDescription) {
		this.referenceDescription = referenceDescription;
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

	public Set<GOterm> getReferenceGoAnnoatations() {
		return referenceGoAnnoatations;
	}

	public void setReferenceGoAnnoatations(Set<GOterm> referenceGoAnnoatations) {
		this.referenceGoAnnoatations = referenceGoAnnoatations;
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
