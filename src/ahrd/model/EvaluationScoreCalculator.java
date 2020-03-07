package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import sun.tools.tree.ThisExpression;

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
	private Set<GOterm> groundTruthGoAnnotations = new HashSet<GOterm>();
	private Fscore simpleGoAnnotationScore;
	private Fscore ancestryGoAnnotationScore;
	private Fscore semSimGoAnnotationScore;
	private Fscore highestPossibleSimpleGoAnnotationScore;
	private Fscore highestPossibleAncestryGoAnnotationScore;
	private Fscore highestPossibleSemSimGoAnnotationScore;
	private Double highestPossibleDescriptionPrecision;
	private Double highestPossibleSimpleGoAnnotationPrecision;
	private Double highestPossibleAncestryGoAnnotationPrecision;
	private Double highestPossibleSemSimGoAnnotationPrecision;
	private Double highestPossibleDescriptionRecall;
	private Double highestPossibleSimpleGoAnnotationRecall;
	private Double highestPossibleAncestryGoAnnotationRecall;
	private Double highestPossibleSemSimGoAnnotationRecall;

	public EvaluationScoreCalculator(Protein protein) {
		super();
		setProtein(protein);
		// In case the ground truth description is empty or has no valid tokens left after the token blacklist is applied,
		// a blast result with an empty description (eg. after the token blacklist is applied) will result in a evaluation score of 1.
		// To make sure the same score is given if no result exists for a database, an empty bestUnchangedBlastResult is set as default for each database.
		// In the vast majority of cases it will be replaced with a proper blast result.
		if (getSettings().getWriteBestBlastHitsToOutput()) {
			for (String blastDbName : getSettings().getBlastDatabases()) {
				BlastResult emptyBlastResult = new BlastResult(blastDbName, "", "");
				emptyBlastResult.setShortAccession("");
				emptyBlastResult.setBitScore(Double.MIN_VALUE);
				getBestUnchangedBlastResults().put(blastDbName, emptyBlastResult);
			}
		}
		// Necessary if the ground truth description is empty or has no valid tokens left after the token blacklist is applied.
		// Also necessary if the ground truth go annotation is empty.
		// In both cases no annotation should result in an evaluation score of 1.
		// To ensure so an empty CompetitorAnnotation is created for each competitor.
		// In the vast majority of cases these CompetitorAnnotations will be replaced by proper annotations read from file.
		if (getSettings().hasCompetitors()) {
			for (String competitor : getSettings().getCompetitorSettings().keySet()) {
				addCompetitorAnnotation(competitor, new CompetitorAnnotation("", ""));
			}
		}
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
		Fscore fBetaScore = new Fscore();
		// Validate Ground Truth:
		if (groundTruthTkns!=null && !groundTruthTkns.isEmpty() && assignedTkns!=null && !assignedTkns.isEmpty()) {
			double tp = truePositives(assignedTkns, groundTruthTkns);
			// Avoid division by zero:
			if (tp > 0.0) {
				fBetaScore.setPrecision(tp / assignedTkns.size());
				fBetaScore.setRecall(tp / groundTruthTkns.size());
			}
		} else {
			if (groundTruthTkns==null || groundTruthTkns.isEmpty()) {
				fBetaScore.setRecall(Double.NaN);
				fBetaScore.setPrecision(Double.NaN);
			}
			if (assignedTkns==null || assignedTkns.isEmpty()) {
				fBetaScore.setPrecision(Double.NaN);
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
			Set<String> hrdEvlTkns;
			if (getProtein().getDescriptionScoreCalculator().getHighestScoringBlastResult() != null) {
				// Generate the set of evaluation-tokens from the actually assigned description.
				// If evaluateValidTokens is set to false: WITHOUT filtering each token with the BLACKLIST.
				getProtein().getDescriptionScoreCalculator().getHighestScoringBlastResult().tokenizeForEvaluation();
				hrdEvlTkns = getProtein().getDescriptionScoreCalculator().getHighestScoringBlastResult().getEvaluationTokens();
			} else {
				// Well, no Description assigned means no tokens:
				hrdEvlTkns = new HashSet<String>();
			}
			// Calculate the Evaluation-Score as the F-Beta-Score (including Precision and Recall):
			setEvalutionScore(fBetaScore(hrdEvlTkns, getGroundTruthDescription().getTokens()));
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
									annot.setSimpleGoAnnotationScore(calcSimpleGoAnnotationScore(this.groundTruthGoAnnotations, annot.getGoAnnotations()));
								if (getSettings().doCalculateAncestryGoF1Scores())
									annot.setAncestryGoAnnotationScore(calcAncestryGoAnnotationScore(this.groundTruthGoAnnotations, annot.getGoAnnotations()));
								if (getSettings().doCalculateSemSimGoF1Scores())
									annot.setSemSimGoAnnotationScore(calcSemSimGoAnnotationScore(this.groundTruthGoAnnotations, annot.getGoAnnotations(), getSettings().doEvaluateSubontologiesSepatatey()));
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
										calcSimpleGoAnnotationScore(this.groundTruthGoAnnotations, cmpt.getGoAnnotations()));
							if (getSettings().doCalculateAncestryGoF1Scores())
								cmpt.setAncestryGoAnnotationScore(
										calcAncestryGoAnnotationScore(this.groundTruthGoAnnotations, cmpt.getGoAnnotations()));
							if (getSettings().doCalculateSemSimGoF1Scores())
								cmpt.setSemSimGoAnnotationScore(
										calcSemSimGoAnnotationScore(this.groundTruthGoAnnotations, cmpt.getGoAnnotations(), getSettings().doEvaluateSubontologiesSepatatey()));
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
				setSimpleGoAnnotationScore(calcSimpleGoAnnotationScore(this.groundTruthGoAnnotations, this.protein.getGoResultsTerms()));
			// Calculation of an F1-score based on ground truth and prediction GO
			// annotations extended to their complete ancestry
			if (getSettings().doCalculateAncestryGoF1Scores())
				setAncestryGoAnnotationScore(calcAncestryGoAnnotationScore(this.groundTruthGoAnnotations, this.protein.getGoResultsTerms()));
			// Calculation of an F1-score based on the semantic similarity
			// (based on term information content) of ground truth and prediction
			// GO annotations.
			if (getSettings().doCalculateSemSimGoF1Scores())
				setSemSimGoAnnotationScore(calcSemSimGoAnnotationScore(this.groundTruthGoAnnotations, this.protein.getGoResultsTerms(), getSettings().doEvaluateSubontologiesSepatatey()));
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
		for (Iterator<GOterm> groundTruthIter = groundTruth.iterator(); groundTruthIter.hasNext();) {
			if (prediction.contains(groundTruthIter.next())) {
				truePositive++;
			}
		}
		f.setRecall((double) truePositive / groundTruth.size());
		if (groundTruth.size() > 0) {
			f.setPrecision((double) truePositive / prediction.size());
		} else {
			f.setPrecision(Double.NaN);
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
		for (Iterator<GOterm> groundTruthIter = groundTruthAncestry.iterator(); groundTruthIter.hasNext();) {
			if (predictionAncestry.contains(groundTruthIter.next())) {
				truePositive++;
			}
		}
		recall = (double) truePositive / groundTruthAncestry.size();
		if (groundTruthAncestry.size() > 0) {
			precision = (double) truePositive / predictionAncestry.size();
		} else {
			precision = Double.NaN;
		}
		f.setPrecision(precision);
		f.setRecall(recall);
		return f;
	}

	/**
	 * Calculates an F-score from recall and precision based on the semantic
	 * similarity of the ground truth and prediction GO term sets.
	 * 
	 * @return Double - F score
	 */
	private Fscore calcSemSimGoAnnotationScore(Set<GOterm> groundTruth, Set<GOterm> prediction) {
		Fscore f = new Fscore();
		if (getSettings().doWriteCumulativeSemSimGoScores()) {
			f = new InfoContentFscore();
		}
		Double recall = 0.0;
		Double precision = 0.0;
		/**
		 * Find the information content of the ground truth: Usually the
		 * information content of the terms themselves, but if they havn't been
		 * used in swissprot their info content is infinite so the highest non
		 * infinite info content from their ancestry is used as fallback.
		 */
		if (groundTruth.size() > 0 && prediction.size() > 0) {
			// Recall
			Double infoContentGroundTruth = 0.0;
			Double commonInfoContentPrediction = 0.0;
			for (GOterm groundTruthTerm : groundTruth) {
				Double maxInfoContentAncestry = 0.0;
				for (GOterm ancestryTerm : groundTruthTerm.getAncestry()) {
					Double infoContent = ancestryTerm.getInformationContent();
					if (Double.isFinite(infoContent) && infoContent > maxInfoContentAncestry) {
						maxInfoContentAncestry = infoContent;
					}
				}
				infoContentGroundTruth += maxInfoContentAncestry;
				Double maxCommonInfoContentPrediction = 0.0;
				for (GOterm predictionTerm : prediction) {
					Double infoContent = maxCommonInformationContent(groundTruthTerm, predictionTerm);
					if (Double.isFinite(infoContent) && infoContent > maxCommonInfoContentPrediction) {
						maxCommonInfoContentPrediction = infoContent;
					}
				}
				commonInfoContentPrediction += maxCommonInfoContentPrediction;
			}
			if (infoContentGroundTruth > 0.0) {
				recall = commonInfoContentPrediction / infoContentGroundTruth;
				if (getSettings().doWriteCumulativeSemSimGoScores()) {
					InfoContentFscore icf = (InfoContentFscore)f;
					icf.setCommonInfoContentPrediction(commonInfoContentPrediction);
					icf.setInfoContentGroundTruth(infoContentGroundTruth);
				}
			} else {
				recall = Double.NaN;
			}
			// Precision
			Double infoContentPrediction = 0.0;
			Double commonInfoContentGroundTruth = 0.0;
			for (GOterm predictionTerm : prediction) {
				Double maxInfoContentAncestry = 0.0;
				for (GOterm ancestryTerm : predictionTerm.getAncestry()) {
					Double infoContent = ancestryTerm.getInformationContent();
					if (Double.isFinite(infoContent) && infoContent > maxInfoContentAncestry) {
						maxInfoContentAncestry = infoContent;
					}
				}
				Double maxCommonInfoContentGroundTruth = -1.0;
				for (GOterm groundTruthTerm : groundTruth) {
					Set<GOterm> commonAncestry = new HashSet<GOterm>(predictionTerm.getAncestry());
					commonAncestry.retainAll(groundTruthTerm.getAncestry());
					/*
					 * If the maxCommonInfoContent remains negative, the commonAncestry is empty.
					 * -> This means the predicted GO term and the ground truth GO term are from different ontologies.
					 * If the maxCommonInfoContent turns out to 0.0, the predicted GO term and the ground truth GO term are from the same ontology but only share the root term in their ancestries.
					 */
					Double maxCommonInfoContent = -1.0;
					for (GOterm ancestryTerm : commonAncestry) {
						Double infoContent = ancestryTerm.getInformationContent();
						if (Double.isFinite(infoContent) && infoContent > maxCommonInfoContent) {
							maxCommonInfoContent = infoContent;
						}
					}
					if (maxCommonInfoContent > maxCommonInfoContentGroundTruth) {
						maxCommonInfoContentGroundTruth = maxCommonInfoContent;
					}
				}
				/*
				 * If the maxCommonInfoContentGroundTruth remains negative, there are no GO terms in the ground truth within the same ontology as the predicted GO term.
				 * -> This means in the ontology of the predicted GO term is no knowledge available for the protein in question.
				 * -> Thus the prediction also can't be proven wrong and should not be penalized.
				 * -> Consequently the predicted GO term's information content will not be added to the information content of the prediction.   
				 * If the maxCommonInfoContentGroundTruth turns out to 0.0, ground truth GO terms were available in the same ontology but the only common ancestor was the root term.
				 * -> Thus information about the function of the proteins was available and the prediction was very wrong.
				 * -> Consequently the predicted GO term's information content will be added to the information content of the prediction while 0 will be added to the information content common with the ground truth.
				 * -> This results in a penalizing of the prediction to the full extent. 
				 */
				if (maxCommonInfoContentGroundTruth >= 0.0) {
					infoContentPrediction += maxInfoContentAncestry;
					commonInfoContentGroundTruth += maxCommonInfoContentGroundTruth;
				}
			}
			if (infoContentPrediction > 0.0) {
				precision = commonInfoContentGroundTruth / infoContentPrediction;
				if (getSettings().doWriteCumulativeSemSimGoScores()) {
					InfoContentFscore icf = (InfoContentFscore)f;
					icf.setCommonInfoContentGroundTruth(commonInfoContentGroundTruth);
					icf.setInfoContentPrediction(infoContentPrediction);
				}
			} else {
				precision = Double.NaN;
			}
		} else {
			if (groundTruth.size() == 0) {
				recall = Double.NaN;
				precision = Double.NaN;
			}
			if (prediction.size() == 0) {
				precision = Double.NaN;
			}
		}
		if (f.getScore() > 1.0) { // Something went very wrong - Should never be happening
			System.out.println("p: " + precision + "\tr: " + recall + "\tf1: " + f);
		}
		f.setPrecision(precision);
		f.setRecall(recall);
		return f;
	}

	private Double maxCommonInformationContent(GOterm firstTerm, GOterm secondTerm) {
		Set<GOterm> commonAncestry = new HashSet<GOterm>(firstTerm.getAncestry());
		commonAncestry.retainAll(secondTerm.getAncestry());
		Double maxCommonInfoContent = 0.0;
		for (GOterm ancestryTerm : commonAncestry) {
			Double infoContent = ancestryTerm.getInformationContent();
			if (Double.isFinite(infoContent) && infoContent > maxCommonInfoContent) {
				maxCommonInfoContent = infoContent;
			}
		}
		return maxCommonInfoContent;
	}
	
	private Fscore calcSemSimGoAnnotationScore(Set<GOterm> groundTruth, Set<GOterm> prediction, boolean perSubOntology) {
		Fscore basicFscore = calcSemSimGoAnnotationScore(groundTruth, prediction);
		if (!perSubOntology) {
			return basicFscore;
		} else {
			GoFscore subOntologyFscore = new GoFscore();
			subOntologyFscore.setAllFscore(basicFscore);
			// BPO
			Set<GOterm> bpoGroundTruth = new HashSet<>();
			for (GOterm term : groundTruth) {
				if (term.getOntology().equals("biological_process")) {
					bpoGroundTruth.add(term);
				}
			}
			Set<GOterm> bpoPrediction = new HashSet<>();
			for (GOterm term : prediction) {
				if (term.getOntology().equals("biological_process")) {
					bpoPrediction.add(term);
				}
			}
			subOntologyFscore.setBpoFscore(calcSemSimGoAnnotationScore(bpoGroundTruth, bpoPrediction));
			// MFO
			Set<GOterm> mfoGroundTruth = new HashSet<>();
			for (GOterm term : groundTruth) {
				if (term.getOntology().equals("molecular_function")) {
					mfoGroundTruth.add(term);
				}
			}
			Set<GOterm> mfoPrediction = new HashSet<>();
			for (GOterm term : prediction) {
				if (term.getOntology().equals("molecular_function")) {
					mfoPrediction.add(term);
				}
			}
			subOntologyFscore.setMfoFscore(calcSemSimGoAnnotationScore(mfoGroundTruth, mfoPrediction));
			// CCO
			Set<GOterm> ccoGroundTruth = new HashSet<>();
			for (GOterm term : groundTruth) {
				if (term.getOntology().equals("cellular_component")) {
					ccoGroundTruth.add(term);
				}
			}
			Set<GOterm> ccoPrediction = new HashSet<>();
			for (GOterm term : prediction) {
				if (term.getOntology().equals("cellular_component")) {
					ccoPrediction.add(term);
				}
			}
			subOntologyFscore.setCcoFscore(calcSemSimGoAnnotationScore(ccoGroundTruth, ccoPrediction));
			
			return subOntologyFscore;
		}
	}

	/**
	 * In order to get more accurate information of how well AHRD performs, we
	 * infer the highest possible score by calculating the evaluation-score for
	 * each BlastResult's Description and remembering the highest achieved
	 * score.
	 */
	public void findBlastResultWithHighestPossibleDescriptionScore() {
		setHighestPossibleDescriptionScore(fBetaScore(new HashSet<String>(), getGroundTruthDescription().getTokens())); // In case the ground truth is empty
		for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
			for (BlastResult cmpt : resultsFromBlastDatabase) {
				// Generate the set of evaluation-tokens for each description, 
				// if evaluateValidTokens is set to false, WITHOUT filtering each token with the BLACKLIST.
				cmpt.tokenizeForEvaluation();
				cmpt.setEvaluationScore(fBetaScore(cmpt.getEvaluationTokens(), getGroundTruthDescription().getTokens()));
				// Find best performing BlastResult-Description:
				if (cmpt.getEvaluationScore().getScore() > getHighestPossibleDescriptionScore().getScore() || getHighestPossibleDescriptionScore().getScore().isNaN()) {
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
				this.setHighestPossibleSimpleGoAnnotationScore(calcSimpleGoAnnotationScore(this.groundTruthGoAnnotations, new HashSet<GOterm>())); // In case ground truth go annotation is empty
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcSimpleGoAnnotationScore(this.groundTruthGoAnnotations, br.getGoAnnotations());
						if (score.getScore() > this.getHighestPossibleSimpleGoAnnotationScore().getScore() || this.getHighestPossibleSimpleGoAnnotationScore().getScore().isNaN())
								this.setHighestPossibleSimpleGoAnnotationScore(score);
						if (getHighestPossibleSimpleGoAnnotationScore().getScore().equals(1.0))
							break;
					}
					if (getHighestPossibleSimpleGoAnnotationScore().getScore().equals(1.0))
						break;
				}
			}
			if (getSettings().doCalculateAncestryGoF1Scores()) {
				this.setHighestPossibleAncestryGoAnnotationScore(this.calcAncestryGoAnnotationScore(this.groundTruthGoAnnotations, new HashSet<GOterm>())); // In case ground truth go annotation is empty
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcAncestryGoAnnotationScore(this.groundTruthGoAnnotations, br.getGoAnnotations());
						if (score.getScore() > this.getHighestPossibleAncestryGoAnnotationScore().getScore() || this.getHighestPossibleAncestryGoAnnotationScore().getScore().isNaN())
								this.setHighestPossibleAncestryGoAnnotationScore(score);
						if (getHighestPossibleAncestryGoAnnotationScore().getScore().equals(1.0))
							break;
					}
					if (getHighestPossibleAncestryGoAnnotationScore().getScore().equals(1.0))
						break;
				}
			}
			if (getSettings().doCalculateSemSimGoF1Scores()) {
				this.setHighestPossibleSemSimGoAnnotationScore(calcSemSimGoAnnotationScore(this.groundTruthGoAnnotations, new HashSet<GOterm>(), getSettings().doEvaluateSubontologiesSepatatey())); // In case ground truth go annotation is empty
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcSemSimGoAnnotationScore(this.groundTruthGoAnnotations, br.getGoAnnotations(), getSettings().doEvaluateSubontologiesSepatatey());
						if (score.getScore() > this.getHighestPossibleSemSimGoAnnotationScore().getScore() || this.getHighestPossibleSemSimGoAnnotationScore().getScore().isNaN())
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

	public void findBlastResultWithHighestPossiblePrecision() {
		// find the blast result with the description that results in the highest possible precision
		if (getGroundTruthDescription().getTokens().isEmpty()) { 
			setHighestPossibleDescriptionPrecision(Double.NaN);
		} else {
			setHighestPossibleDescriptionPrecision(0.0);
			for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
				for (BlastResult br : resultsFromBlastDatabase) {
					// Generate the set of evaluation-tokens for each description, 
					// if evaluateValidTokens is set to false, WITHOUT filtering each token with the BLACKLIST.
					br.tokenizeForEvaluation();
					br.setEvaluationScore(fBetaScore(br.getEvaluationTokens(), getGroundTruthDescription().getTokens()));
					// Find best performing BlastResult-Description:
					if (br.getEvaluationScore().getPrecision() > getHighestPossibleDescriptionPrecision()) {
						setHighestPossibleDescriptionPrecision(br.getEvaluationScore().getPrecision());
					}
					if (getHighestPossibleDescriptionPrecision().equals(1.0))
						break;
				}
			}
		}
		// find the blast result with the go annotations that results in the highest possible precision
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			if (getSettings().doCalculateSimpleGoF1Scores()) {
				if(this.groundTruthGoAnnotations.isEmpty()) {
					this.setHighestPossibleSimpleGoAnnotationPrecision(Double.NaN);
				} else {
					this.setHighestPossibleSimpleGoAnnotationPrecision(0.0);
					for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
						for (BlastResult br : resultsFromBlastDatabase) {
							Fscore score = calcSimpleGoAnnotationScore(this.groundTruthGoAnnotations, br.getGoAnnotations());
							if (score.getPrecision() > this.getHighestPossibleSimpleGoAnnotationPrecision())
									this.setHighestPossibleSimpleGoAnnotationPrecision(score.getPrecision());
							if (getHighestPossibleSimpleGoAnnotationPrecision().equals(1.0))
								break;
						}
						if (getHighestPossibleSimpleGoAnnotationPrecision().equals(1.0))
							break;
					}
				}
			}
			if (getSettings().doCalculateAncestryGoF1Scores()) {
				if(this.groundTruthGoAnnotations.isEmpty()) {
					this.setHighestPossibleAncestryGoAnnotationPrecision(Double.NaN);
				} else {
					this.setHighestPossibleAncestryGoAnnotationPrecision(0.0);
					for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
						for (BlastResult br : resultsFromBlastDatabase) {
							Fscore score = calcAncestryGoAnnotationScore(this.groundTruthGoAnnotations, br.getGoAnnotations());
							if (score.getPrecision() > this.getHighestPossibleAncestryGoAnnotationPrecision())
									this.setHighestPossibleAncestryGoAnnotationPrecision(score.getPrecision());
							if (getHighestPossibleAncestryGoAnnotationPrecision().equals(1.0))
								break;
						}
						if (getHighestPossibleAncestryGoAnnotationPrecision().equals(1.0))
							break;
					}
				}
			}
			if (getSettings().doCalculateSemSimGoF1Scores()) {
				if (this.groundTruthGoAnnotations.isEmpty()) {
					this.setHighestPossibleSemSimGoAnnotationPrecision(Double.NaN);
				} else {
					this.setHighestPossibleSemSimGoAnnotationPrecision(0.0);
					for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
						for (BlastResult br : resultsFromBlastDatabase) {
							Fscore score = calcSemSimGoAnnotationScore(this.groundTruthGoAnnotations, br.getGoAnnotations());
							if (score.getPrecision() > this.getHighestPossibleSemSimGoAnnotationPrecision())
									this.setHighestPossibleSemSimGoAnnotationPrecision(score.getPrecision());
							if (getHighestPossibleSemSimGoAnnotationPrecision().equals(1.0))
								break;
						}
						if (getHighestPossibleSemSimGoAnnotationPrecision().equals(1.0))
							break;
					}
				}
			}
		}
	}

	public void findBlastResultWithHighestPossibleRecall() {
		// find the blast result with the description that results in the highest possible recall 
		setHighestPossibleDescriptionRecall(fBetaScore(new HashSet<String>(), getGroundTruthDescription().getTokens()).getRecall()); // NaN in case ground truth description is empty, initialize with 0 otherwise
		for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
			for (BlastResult br : resultsFromBlastDatabase) {
				// Generate the set of evaluation-tokens for each description, 
				// if evaluateValidTokens is set to false, WITHOUT filtering each token with the BLACKLIST.
				br.tokenizeForEvaluation();
				br.setEvaluationScore(fBetaScore(br.getEvaluationTokens(), getGroundTruthDescription().getTokens()));
				// Find best performing BlastResult-Description:
				if (br.getEvaluationScore().getRecall() > getHighestPossibleDescriptionRecall()) {
					setHighestPossibleDescriptionRecall(br.getEvaluationScore().getRecall());
				}
				if (getHighestPossibleDescriptionRecall().equals(1.0))
					break;
			}
		}
		// find the blast result with the go annotations that results in the highest possible recall
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			if (getSettings().doCalculateSimpleGoF1Scores()) {
				this.setHighestPossibleSimpleGoAnnotationRecall(calcSimpleGoAnnotationScore(this.groundTruthGoAnnotations, new HashSet<GOterm>()).getRecall()); // NaN in case ground truth go annotation is empty, initialize with 0 otherwise
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcSimpleGoAnnotationScore(this.groundTruthGoAnnotations, br.getGoAnnotations());
						if (score.getRecall() > this.getHighestPossibleSimpleGoAnnotationRecall())
								this.setHighestPossibleSimpleGoAnnotationRecall(score.getRecall());
						if (getHighestPossibleSimpleGoAnnotationRecall().equals(1.0))
							break;
					}
					if (getHighestPossibleSimpleGoAnnotationRecall().equals(1.0))
						break;
				}
			}
			if (getSettings().doCalculateAncestryGoF1Scores()) {
				this.setHighestPossibleAncestryGoAnnotationRecall(calcAncestryGoAnnotationScore(this.groundTruthGoAnnotations, new HashSet<GOterm>()).getRecall()); // NaN in case ground truth go annotation is empty, initialize with 0 otherwise
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcAncestryGoAnnotationScore(this.groundTruthGoAnnotations, br.getGoAnnotations());
						if (score.getRecall() > this.getHighestPossibleAncestryGoAnnotationRecall())
								this.setHighestPossibleAncestryGoAnnotationRecall(score.getRecall());
						if (getHighestPossibleAncestryGoAnnotationRecall().equals(1.0))
							break;
					}
					if (getHighestPossibleAncestryGoAnnotationRecall().equals(1.0))
						break;
				}
			}
			if (getSettings().doCalculateSemSimGoF1Scores()) {
				this.setHighestPossibleSemSimGoAnnotationRecall(calcSemSimGoAnnotationScore(this.groundTruthGoAnnotations, new HashSet<GOterm>()).getRecall()); // NaN in case ground truth go annotation is empty, initialize with 0 otherwise  
				for (List<BlastResult> resultsFromBlastDatabase : getProtein().getBlastResults().values()) {
					for (BlastResult br : resultsFromBlastDatabase) {
						Fscore score = calcSemSimGoAnnotationScore(this.groundTruthGoAnnotations, br.getGoAnnotations());
						if (score.getRecall() > this.getHighestPossibleSemSimGoAnnotationRecall())
								this.setHighestPossibleSemSimGoAnnotationRecall(score.getRecall());
						if (getHighestPossibleSemSimGoAnnotationRecall().equals(1.0))
							break;
					}
					if (getHighestPossibleSemSimGoAnnotationRecall().equals(1.0))
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

	public Set<GOterm> getGroundTruthGoAnnotations() {
		return groundTruthGoAnnotations;
	}

	public void setGroundTruthGoAnnotations(Set<GOterm> groundTruthGoAnnotations) {
		this.groundTruthGoAnnotations = groundTruthGoAnnotations;
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

	public Double getHighestPossibleDescriptionPrecision() {
		return highestPossibleDescriptionPrecision;
	}

	public void setHighestPossibleDescriptionPrecision(Double highestPossibleDescriptionPrecision) {
		this.highestPossibleDescriptionPrecision = highestPossibleDescriptionPrecision;
	}

	public Double getHighestPossibleSimpleGoAnnotationPrecision() {
		return highestPossibleSimpleGoAnnotationPrecision;
	}

	public void setHighestPossibleSimpleGoAnnotationPrecision(Double highestPossibleSimpleGoAnnotationPrecision) {
		this.highestPossibleSimpleGoAnnotationPrecision = highestPossibleSimpleGoAnnotationPrecision;
	}

	public Double getHighestPossibleAncestryGoAnnotationPrecision() {
		return highestPossibleAncestryGoAnnotationPrecision;
	}

	public void setHighestPossibleAncestryGoAnnotationPrecision(Double highestPossibleAncestryGoAnnotationPrecision) {
		this.highestPossibleAncestryGoAnnotationPrecision = highestPossibleAncestryGoAnnotationPrecision;
	}

	public Double getHighestPossibleSemSimGoAnnotationPrecision() {
		return highestPossibleSemSimGoAnnotationPrecision;
	}

	public void setHighestPossibleSemSimGoAnnotationPrecision(Double highestPossibleSemSimGoAnnotationPrecision) {
		this.highestPossibleSemSimGoAnnotationPrecision = highestPossibleSemSimGoAnnotationPrecision;
	}

	public Double getHighestPossibleDescriptionRecall() {
		return highestPossibleDescriptionRecall;
	}

	public void setHighestPossibleDescriptionRecall(Double highestPossibleDescriptionRecall) {
		this.highestPossibleDescriptionRecall = highestPossibleDescriptionRecall;
	}

	public Double getHighestPossibleSimpleGoAnnotationRecall() {
		return highestPossibleSimpleGoAnnotationRecall;
	}

	public void setHighestPossibleSimpleGoAnnotationRecall(Double highestPossibleSimpleGoAnnotationRecall) {
		this.highestPossibleSimpleGoAnnotationRecall = highestPossibleSimpleGoAnnotationRecall;
	}

	public Double getHighestPossibleAncestryGoAnnotationRecall() {
		return highestPossibleAncestryGoAnnotationRecall;
	}

	public void setHighestPossibleAncestryGoAnnotationRecall(Double highestPossibleAncestryGoAnnotationRecall) {
		this.highestPossibleAncestryGoAnnotationRecall = highestPossibleAncestryGoAnnotationRecall;
	}

	public Double getHighestPossibleSemSimGoAnnotationRecall() {
		return highestPossibleSemSimGoAnnotationRecall;
	}

	public void setHighestPossibleSemSimGoAnnotationRecall(Double highestPossibleSemSimGoAnnotationRecall) {
		this.highestPossibleSemSimGoAnnotationRecall = highestPossibleSemSimGoAnnotationRecall;
	}

}
