package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.semanticweb.owlapi.model.OWLOntologyCreationException;

import ahrd.exception.MissingAccessionException;
import ahrd.model.BlastResult;
import ahrd.model.CompetitorAnnotation;
import ahrd.model.GOdatabase;
import ahrd.model.GOterm;
import ahrd.model.Protein;
import ahrd.model.GroundTruthDescription;
import ahrd.view.EvaluatorOutputWriter;
import ahrd.view.TsvOutputWriter;

public class Evaluator extends AHRD {
	
	public Evaluator(String pathToInputYml) throws IOException {
		super(pathToInputYml);
	}

	public void setupGroundTruthDescriptions() throws IOException, MissingAccessionException {
		List<String> fastaEntries = Protein.splitFasta(getSettings().getGroundTruthFasta());
		for (String fastaEntry : fastaEntries) {
			if (fastaEntry != null && !fastaEntry.trim().equals("")) {
				GroundTruthDescription rd = GroundTruthDescription.constructFromFastaEntry(fastaEntry.trim());
				Protein p = getProteins().get(rd.getAccession());
				if (p == null)
					throw new MissingAccessionException(
							"Could not find Protein for Accession '" + rd.getAccession() + "'");
				p.getEvaluationScoreCalculator().setGroundTruthDescription(rd);
			}
		}
	}

	public void setupGoAnnotationEvaluation() throws OWLOntologyCreationException, IOException, MissingAccessionException {
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
			// Load a Map of all GO terms
			if (goDB == null) {
				goDB = new GOdatabase().getMap();
			}
			// Load ground truth GO annotations
			for (String groundTruthGoAnnotationFileEntryLine : getSettings().getGroundTruthGoAnnotationsFromFile()) {
				String[] groundTruthGoAnnotationFileEntry = groundTruthGoAnnotationFileEntryLine.split("\t");
				String protAcc = groundTruthGoAnnotationFileEntry[0].trim();
				String termAcc = groundTruthGoAnnotationFileEntry[1].trim();
				Protein p = getProteins().get(protAcc);
				if (p == null) {
					p = getShortAccsProteins().get(protAcc); // Otherwise try if the ground truth uses short protein accessions
					if (p == null) {
						throw new MissingAccessionException("Could not find protein for accession '" + protAcc + "'");
					}
				}
				GOterm term = goDB.get(termAcc);
				if (term == null) {
					throw new MissingAccessionException("Could not find GO term for accession '" + termAcc + "'");
				}
				p.getEvaluationScoreCalculator().getGroundTruthGoAnnotations().add(term);
			}
			// Add GOterm objects to predicted annotations
			goAnnotsStringToObject();
			// Annotate the best blast results with GOterm objects
			if (getSettings().getWriteBestBlastHitsToOutput()) {
				for (Iterator<Protein> protIter = getProteins().values().iterator(); protIter.hasNext();){
					Map<String, BlastResult> bestBlastResult = protIter.next().getEvaluationScoreCalculator().getBestUnchangedBlastResults();
					if (bestBlastResult != null) {
						for (String blastDb : bestBlastResult.keySet()) {
						String bestBlastResultShortAccession = bestBlastResult.get(blastDb).getShortAccession();
							if (getGoAnnotationReference().containsKey(bestBlastResultShortAccession)) {
								for (String termAcc : getGoAnnotationReference().get(bestBlastResultShortAccession)) {
									GOterm term = goDB.get(termAcc);
									if (term == null) {
										throw new MissingAccessionException("Could not find GO term for accession '" + termAcc + "'");
									}
									bestBlastResult.get(blastDb).getGoAnnotations().add(term);
								}
							}
						}
					}
				}
			}
			// If highest possible GO-annotation scores are requested all BlastResults have to be annotated first
			if (getSettings().doFindHighestPossibleGoScore() || getSettings().doFindHighestPossiblePrecision() || getSettings().doFindHighestPossibleRecall()) {
				for (Protein prot : getProteins().values()) {
					for (List<BlastResult> blastDbResults : prot.getBlastResults().values()) {
						for (BlastResult br : blastDbResults) {
							String blastResultShortAccession = br.getShortAccession();
							if (getGoAnnotationReference().containsKey(blastResultShortAccession)) {
								for (String termAcc : getGoAnnotationReference().get(blastResultShortAccession)) {
									GOterm term = goDB.get(termAcc);
									if (term == null) {
										throw new MissingAccessionException("Could not find GO term for accession '" + termAcc + "'");
									}
									br.getGoAnnotations().add(term);
								}
							}
						}
					}
				}
			}
		}
	}

	public void goAnnotsStringToObject() throws MissingAccessionException, OWLOntologyCreationException, IOException {
		// Load a Map of all GO terms
		if (goDB == null)
			goDB = new GOdatabase().getMap();
		// Add GOterm objects to predicted annotations
		for (Iterator<Protein> protIter = getProteins().values().iterator(); protIter.hasNext();) {
			Protein prot = protIter.next();
			prot.setGoResultsTerms(new HashSet<GOterm>());
			for (String termAcc : prot.getGoResults()) {
				GOterm term = goDB.get(termAcc);
				if (term == null) {
					throw new MissingAccessionException("Could not find GO term for accession '" + termAcc + "'");
				}
				prot.getGoResultsTerms().add(term);
			}
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Usage:\njava -Xmx2g -cp ahrd.jar ahrd.controller.Evaluator input.yml\n");

		try {
			Evaluator evaluator = new Evaluator(args[0]);
			evaluator.setup(false); // false -> Don't log memory and time-usages
			// After the setup the unique short accessions are no longer needed:
			evaluator.setUniqueBlastResultShortAccessions(null);
			evaluator.setupGroundTruthDescriptions();
			// Iterate over all Proteins and assign the best scoring Human
			// Readable Description
			evaluator.assignHumanReadableDescriptions();
			// Load a Map of all GO terms
			// Load ground truth GO annotations
			// Add GOterm objects to predicted annotations
			// Annotate the best blast results with GOterm objects
			evaluator.setupGoAnnotationEvaluation();
			// Sets up description and GOterm annotations of competitors in the EvaluationScoreCalculators of each Protein   
			evaluator.setupCompetitors();
			// Evaluate AHRD's and Competitors Performance for each Protein:
			evaluator.calculateEvaluationScores();
			// If requested, calculate the highest possibly achievable
			// evaluation score:
			if (getSettings().doFindHighestPossibleEvaluationScore())
				evaluator.findHighestPossibleEvaluationScores();
			// If requested, calculate the highest possibly achievable score for go annotations:
			if (getSettings().doFindHighestPossibleGoScore())
				evaluator.findHighestPossibleGoScores();
			// If requested, calculate the highest possibly achievable precision for descriptions and go annotations:
			if (getSettings().doFindHighestPossiblePrecision())
				evaluator.findHighestPossiblePrecision();
			// If requested, calculate the highest possibly achievable recall for descriptions and go annotations:
			if (getSettings().doFindHighestPossibleRecall())
				evaluator.findHighestPossibleRecall();
			// Generate Output:
			TsvOutputWriter ow = new EvaluatorOutputWriter(evaluator.getProteins().values());
			ow.writeOutput();
			System.out.println("Written output into:\n" + getSettings().getPathToOutput());
		} catch (Exception e) {
			System.err.println("We are sorry, an unexpected ERROR occurred:");
			e.printStackTrace(System.err);
		}

	}

	/**
	 * Calculates the evaluation-score for the description assigned by AHRD and
	 * all found descriptions coming from competitive methods. As for now those
	 * methods are the best Blast-Hits of each searched Blast-Database. @NOTE:
	 * Best Blast-Hits have the highest Bit-Scores. <br />
	 * After having evaluated the objective function for each Protein,
	 * calculates the average-evaluation-score and sets it in the current
	 * Settings.
	 */
	public void calculateEvaluationScores() {
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator().assignEvaluationScores();
		}
	}

	/**
	 * Calculates the evaluation-score as defined in function
	 * calculateEvaluationScores() for each Protein's Blast-Results'
	 * Description-Lines and stores the maximum possible score each Protein's
	 * EvaluationScoreCalculator. This is used to get more accurate information
	 * of how well AHRD performs.
	 */
	public void findHighestPossibleEvaluationScores() {
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator().findBlastResultWithHighestPossibleDescriptionScore();
		}
	}
	
	/**
	 * Calculates the go-annotation-scores requested (simple, ancestry, semsim)
	 * for each Protein's Blast-Results' GOterm annotations and stores the 
	 * maximum possible score in each Protein's EvaluationScoreCalculator. 
	 * This is used to get more accurate information of how well AHRD performs.
	 */
	public void findHighestPossibleGoScores() {
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator().findHighestPossibleGoScore();
		}
	}
	
	private void findHighestPossiblePrecision() {
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator().findBlastResultWithHighestPossiblePrecision();
		}		
	}

	private void findHighestPossibleRecall() {
		for (Protein prot : getProteins().values()) {
			prot.getEvaluationScoreCalculator().findBlastResultWithHighestPossibleRecall();
		}		
	}
	
	/**
	 * Annotation files must be tab separated and w/o column headers. 
	 * @throws IOException 
	 * @throws MissingAccessionException 
	 */
	public void setupCompetitors() throws IOException, MissingAccessionException {
		if(getSettings().hasCompetitors()){
			for (String competitor : getSettings().getCompetitorSettings().keySet()) {
				Map<String, CompetitorAnnotation> annots = new HashMap<String, CompetitorAnnotation>();
				// Parse description lines and build a map of protein accessions to description annotations
				for (String competitorDescriptionLine : getSettings().getCompetitorDescriptions(competitor)) {
					String[] values = competitorDescriptionLine.split("\t");
					String accession = values[0].trim();
					CompetitorAnnotation annot = new CompetitorAnnotation(accession, values[1].trim());
					annots.put(accession, annot);
				}
				// GO annotations
				if(getSettings().hasGeneOntologyAnnotations() && getSettings().hasGroundTruthGoAnnotations()) {
					// Parse GO annotation lines and add GO annotations to competitor annotations
					for (String competitorGoAnnotationLine : getSettings().getCompetitorGOAnnotations(competitor)) {
						String[] values = competitorGoAnnotationLine.split("\t");
						String accession = values[0].trim();
						CompetitorAnnotation annot = annots.get(accession);
						if (annot == null) {
							annot = new CompetitorAnnotation(accession, "");
							annots.put(accession, annot);
						}
						String termAcc = values[1].trim();
						GOterm term = this.goDB.get(termAcc);
						if (term == null)
							throw new MissingAccessionException("Could not find '" + termAcc + "' in GOterm database for " + competitor + " Annotation of protein '" + accession + "'");
						annot.getGoAnnotations().add(term);
					}
				}
				// Add competitor annotations to the evaluation score calculators of the appropriate proteins
				for (Map.Entry<String, CompetitorAnnotation> annot : annots.entrySet()) {
					Protein p = getProteins().get(annot.getKey());
					if (p == null)
						throw new MissingAccessionException(
								"Could not find Protein for Accession '" + annot.getKey() + "'");
					p.getEvaluationScoreCalculator().addCompetitorAnnotation(competitor, annot.getValue());
				}	
			}
		}
	}

}