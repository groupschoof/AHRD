package ahrd.controller;

import static ahrd.controller.Settings.getSettings;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;

import ahrd.exception.MissingAccessionException;
import ahrd.model.Blast2GoAnnot;
import ahrd.model.GOdatabase;
import ahrd.model.GOterm;
import ahrd.model.Protein;
import ahrd.model.ReferenceDescription;
import ahrd.view.OutputWriter;

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
		if (getSettings().getPathToBlast2GoAnnotations() != null
				&& !getSettings().getPathToBlast2GoAnnotations().equals("")) {
			for (String blast2GoResultEntry : getSettings().getBlast2GoAnnotations()) {
				Blast2GoAnnot b2ga = Blast2GoAnnot.fromBlast2GoEntry(blast2GoResultEntry);
				Protein p = getProteins().get(b2ga.getAccession());
				if (p == null)
					throw new MissingAccessionException(
							"Could not find Protein for Accession '" + b2ga.getAccession() + "'");
				p.getEvaluationScoreCalculator().addBlast2GoAnnot(b2ga);
			}
		}
	}
	
	public void setupGoAnnotationEvaluation() throws FileNotFoundException, IOException , MissingAccessionException {
		if (getSettings().hasGeneOntologyAnnotations()
				&& getSettings().getPathToReferenceGoAnnotations() != null
				&& new File(getSettings().getPathToReferenceGoAnnotations()).exists()) {
			// Load a Map of all GO terms
			goDB = new GOdatabase().getMap();
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
				p.getEvaluationScoreCalculator().addReferenceGoAnnotation(term);
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
			// Blast2GO is another competitor in the field of annotation of
			// predicted Proteins. AHRD might be compared with B2Gs performance:
			evaluator.setupBlast2GoAnnots();
			// Iterate over all Proteins and assign the best scoring Human
			// Readable Description
			evaluator.assignHumanReadableDescriptions();
			// Load a Map of all GO terms
			// Load reference GO annotations
			// Add GOterm objects to predicted annotations
			evaluator.setupGoAnnotationEvaluation();
			// Evaluate AHRD's performance for each Protein:
			evaluator.calculateEvaluationScores();
			// If requested, calculate the highest possibly achievable
			// evaluation score:
			if (getSettings().doFindHighestPossibleEvaluationScore())
				evaluator.findHighestPossibleEvaluationScores();
			// Generate Output:
			OutputWriter ow = new OutputWriter(evaluator.getProteins().values());
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