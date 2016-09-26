package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ahrd.exception.MissingAccessionException;
import ahrd.model.Blast2GoAnnot;
import ahrd.model.BlastResult;
import ahrd.model.GOdatabase;
import ahrd.model.GOterm;
import ahrd.model.InterproResult;
import ahrd.model.Protein;
import ahrd.model.ReferenceDescription;
import ahrd.view.EvaluatorOutputWriter;
import ahrd.view.TsvOutputWriter;

public class Evaluator extends AHRD {
	
	public Evaluator(String pathToInputYml) throws IOException {
		super(pathToInputYml);
	}

	public void setupReferenceDescriptions() throws IOException, MissingAccessionException {
		List<String> fastaEntries = Protein.splitFasta(getSettings().getReferencesFasta());
		for (String fastaEntry : fastaEntries) {
			if (fastaEntry != null && !fastaEntry.trim().equals("")) {
				ReferenceDescription rd = ReferenceDescription.constructFromFastaEntry(fastaEntry.trim());
				Protein p = getProteins().get(rd.getAccession());
				if (p == null)
					throw new MissingAccessionException(
							"Could not find Protein for Accession '" + rd.getAccession() + "'");
				p.getEvaluationScoreCalculator().setReferenceDescription(rd);
			}
		}
	}

	public void setupBlast2GoAnnots() throws IOException, MissingAccessionException {
		if (getSettings().hasBlast2GoAnnotations()) {
			if(!getSettings().hasGeneOntologyAnnotations() || !getSettings().hasReferenceGoAnnotations()) {
				for (String blast2GoResultEntry : getSettings().getBlast2GoAnnotations()) {
					Blast2GoAnnot b2ga = Blast2GoAnnot.fromBlast2GoEntry(blast2GoResultEntry);
					if (b2ga != null) {
						Protein p = getProteins().get(b2ga.getAccession());
						if (p == null)
							throw new MissingAccessionException(
									"Could not find Protein for Accession '" + b2ga.getAccession() + "'");
						p.getEvaluationScoreCalculator().addBlast2GoAnnot(b2ga);
					}
				}
			} else {
				Map<String, Blast2GoAnnot> annots = new HashMap<String, Blast2GoAnnot>();
				// Parse description lines and build a map of protein accessions and blast2go annotations
				for (String blast2GoResultEntry : getSettings().getBlast2GoAnnotations()) {
					Pattern p = Settings.getBlast2GoAnnotationFileDesclineRegex();
					Matcher m = p.matcher(blast2GoResultEntry);
					if (m.find()) {
						Blast2GoAnnot annot = new Blast2GoAnnot(m.group("shortAccession"), m.group("description")); 
						annot.getGoAnnotations().add(this.goDB.get(m.group("goTerm")));
						annots.put(m.group("shortAccession"), annot);
					}
				}
				// Parse GO annotation lines and add go annotations to blast2go annotations
				for (String blast2GoResultEntry : getSettings().getBlast2GoAnnotations()) {
					Pattern p = Settings.getBlast2GoAnnotationFileAnnotlineRegex();
					Matcher m = p.matcher(blast2GoResultEntry);
					if (m.find()) {
						Blast2GoAnnot annot = annots.get(m.group("shortAccession"));
						if (annot == null)
							throw new MissingAccessionException(
									"Could not find B2GO Annotation for Accession '" + m.group("shortAccession") + "'");
						annot.getGoAnnotations().add(this.goDB.get(m.group("goTerm")));
					}
				}
				// Add blast2go annotations to the evaluation score calculators of the appropriate proteins
				for (Map.Entry<String, Blast2GoAnnot> b2gaMapEntry : annots.entrySet()) {
					Protein p = getProteins().get(b2gaMapEntry.getKey());
					if (p == null)
						throw new MissingAccessionException(
								"Could not find Protein for Accession '" + b2gaMapEntry.getKey() + "'");
					p.getEvaluationScoreCalculator().addBlast2GoAnnot(b2gaMapEntry.getValue());
				}
			}
		}
	}

	public void setupGoAnnotationEvaluation() throws FileNotFoundException, IOException, MissingAccessionException {
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasReferenceGoAnnotations()) {
			// Load a Map of all GO terms
			if (goDB == null) {
				goDB = new GOdatabase().getMap();
			}
			// Load reference GO annotations
			for (String referenceGoAnnotationFileEntryLine : getSettings().getReferenceGoAnnotationsFromFile()) {
				String[] referenceGoAnnotationFileEntry = referenceGoAnnotationFileEntryLine.split("\t");
				String protAcc = referenceGoAnnotationFileEntry[0].trim();
				String termAcc = referenceGoAnnotationFileEntry[1].trim();
				Protein p = getProteins().get(protAcc);
				if (p == null) {
					throw new MissingAccessionException("Could not find protein for accession '" + protAcc + "'");
				}
				GOterm term = goDB.get(termAcc);
				if (term == null) {
					throw new MissingAccessionException("Could not find GO term for accession '" + termAcc + "'");
				}
				p.getEvaluationScoreCalculator().getReferenceGoAnnoatations().add(term);
			}
			// Add GOterm objects to predicted annotations
			for (Iterator<Protein> protIter = getProteins().values().iterator(); protIter.hasNext();) {
				Protein prot = protIter.next();
				for (String termAcc : prot.getGoResults()) {
					GOterm term = goDB.get(termAcc);
					if (term == null) {
						throw new MissingAccessionException("Could not find GO term for accession '" + termAcc + "'");
					}
					prot.getGoResultsTerms().add(term);
				}
			}
			// Annotate the best blast results with GOterm objects
			if (getSettings().getWriteBestBlastHitsToOutput()) {
				for (Iterator<Protein> protIter = getProteins().values().iterator(); protIter.hasNext();){
					Map<String, BlastResult> bestBlastResult = protIter.next().getEvaluationScoreCalculator().getUnchangedBlastResults();
					if (bestBlastResult != null) {
						for (String blastDb : bestBlastResult.keySet()) {
						String bestBlastResultShortAccession = bestBlastResult.get(blastDb).getShortAccession();
							if (getDatabaseGoAnnotations().containsKey(bestBlastResultShortAccession)) {
								for (String termAcc : getDatabaseGoAnnotations().get(bestBlastResultShortAccession)) {
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
			evaluator.setupReferenceDescriptions();
			// Iterate over all Proteins and assign the best scoring Human
			// Readable Description
			evaluator.assignHumanReadableDescriptions();
			// Load a Map of all GO terms
			// Load reference GO annotations
			// Add GOterm objects to predicted annotations
			// Annotate the best blast results with GOterm objects
			evaluator.setupGoAnnotationEvaluation();
			// Blast2GO is another competitor in the field of annotation of
			// predicted Proteins. AHRD might be compared with B2Gs performance:
			evaluator.setupBlast2GoAnnots();
			// Evaluate AHRD's performance for each Protein:
			evaluator.calculateEvaluationScores();
			// If requested, calculate the highest possibly achievable
			// evaluation score:
			if (getSettings().doFindHighestPossibleEvaluationScore())
				evaluator.findHighestPossibleEvaluationScores();
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
			prot.getEvaluationScoreCalculator().assignEvlScrsToCompetitors();
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
			prot.getEvaluationScoreCalculator().findHighestPossibleEvaluationScore();
		}
	}

}