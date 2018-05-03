package ahrd.controller;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Settings.setSettings;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.xml.sax.SAXException;

import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingInterproResultException;
import ahrd.exception.MissingProteinException;
import ahrd.model.BlastResult;
import ahrd.model.GOdatabase;
import ahrd.model.GOterm;
import ahrd.model.InterproResult;
import ahrd.model.Protein;
import ahrd.model.TokenScoreCalculator;
import ahrd.model.ReferenceGoAnnotation;
import ahrd.view.FastaOutputWriter;
import ahrd.view.OutputWriter;
import ahrd.view.TsvOutputWriter;
import nu.xom.ParsingException;

public class AHRD {

	public static final String VERSION = "3.11";

	private Map<String, Protein> proteins;
	private Map<String, Double> descriptionScoreBitScoreWeights = new HashMap<String, Double>();
	private Map<String, Set<ReferenceGoAnnotation>> goAnnotationReference;
	private Set<String> uniqueBlastResultShortAccessions;
	private long timestamp;
	private long memorystamp;
	protected Map<String, GOterm> goDB;
	protected Set<GOterm> goCentricTerms = new HashSet<GOterm>();

	protected long takeTime() {
		// Measure time:
		long now = (new Date()).getTime();
		long measuredSeconds = (now - this.timestamp) / 1000;
		this.timestamp = now;
		return measuredSeconds;
	}

	protected long takeMemoryUsage() {
		// Measure Memory-Usage:
		Runtime rt = Runtime.getRuntime();
		this.memorystamp = (rt.totalMemory() - rt.freeMemory()) / (1024 * 1024);
		return this.memorystamp;
	}

	public static void main(String[] args) {
		System.out.println("Usage:\njava -Xmx2g -jar ahrd.jar input.yml\n");

		try {
			AHRD ahrd = new AHRD(args[0]);
			// Load and parse all inputs
			ahrd.setup(true);
			// After the setup the unique short accessions are no longer needed:
			ahrd.setUniqueBlastResultShortAccessions(null);

			// Iterate over all Proteins and assign the best scoring Human
			// Readable Description
			ahrd.assignHumanReadableDescriptions();
			// If requested iterate over all Proteins and assign the best scoring Gene Ontology terms
			if (getSettings().hasGeneOntologyAnnotations()) {
				ahrd.assignGeneOntologyTerms();
			}
			// Log
			System.out.println("...assigned highestest scoring human readable descriptions in " + ahrd.takeTime()
					+ "sec, currently occupying " + ahrd.takeMemoryUsage() + " MB");
			// If requested assign GO slim terms in addition to detailed GO term annotation
			ahrd.annotateWithGoSlim();
			// If requested term centric annotation is performed as well
			ahrd.termCentricAnnotation();
			// Write result to output-file:
			System.out.println("Writing output to '" + getSettings().getPathToOutput() + "'.");
			OutputWriter ow = initializeOutputWriter(ahrd.getProteins().values());
			ow.writeOutput();
			// Log
			System.out.println("Wrote output in " + ahrd.takeTime() + "sec, currently occupying "
					+ ahrd.takeMemoryUsage() + " MB");

			System.out.println("\n\nDONE");
		} catch (Exception e) {
			System.err.println("We are sorry, an un-expected ERROR occurred:");
			e.printStackTrace(System.err);
		}
	}

	public static OutputWriter initializeOutputWriter(Collection<Protein> proteins) {
		OutputWriter ow = null;
		if (getSettings().doOutputFasta())
			ow = new FastaOutputWriter(proteins);
		else
			ow = new TsvOutputWriter(proteins);
		return ow;
	}

	/**
	 * Constructor initializes this run's settings as a thread-local variable.
	 * Also conditionally initializes fields
	 * <code>uniqueBlastResultShortAccessions</code> and
	 * <code>goAnnotationReference</code> required only if AHRD is requested to
	 * generate Gene Ontology term annotations.
	 * 
	 * @param pathToYmlInput
	 * @throws IOException
	 */
	public AHRD(String pathToYmlInput) throws IOException {
		super();
		setSettings(new Settings(pathToYmlInput));
		// The following fields are only used if AHRD is requested to generate
		// Gene Ontology term annotations:
		if (getSettings().hasGeneOntologyAnnotations()) {
			this.setUniqueBlastResultShortAccessions(new HashSet<String>());
			this.setGoAnnotationReference(new HashMap<String, Set<ReferenceGoAnnotation>>());
		}
	}

	public void initializeProteins() throws IOException, MissingAccessionException {
		setProteins(Protein.initializeProteins(getSettings().getProteinsFasta()));
	}

	public void parseBlastResults() throws IOException, MissingProteinException, MissingAccessionException, SAXException {
		for (String blastDatabase : getSettings().getBlastDatabases()) {
			BlastResult.readBlastResults(getProteins(), blastDatabase, getUniqueBlastResultShortAccessions());
		}
	}

	public void parseInterproResult() throws IOException {
		if (getSettings().hasInterproAnnotations()) {
			Set<String> missingProteinAccessions = new HashSet<String>();
			try {
				InterproResult.parseInterproResult(proteins);
			} catch (MissingProteinException mpe) {
				missingProteinAccessions.add(mpe.getMessage());
			}
			if (missingProteinAccessions.size() > 0) {
				System.err
						.println("WARNING - The following Gene-Accessions were referenced in the Interpro-Result-File,"
								+ " but could not be found in the Memory-Protein-Database:\n"
								+ missingProteinAccessions);
			}
		}
	}

	/**
	 * Method finds GO term annotations for Proteins in the searched Blast
	 * databases and stores them in a Map.
	 * 
	 * @throws IOException
	 */
	public void setUpGoAnnotationReference() throws IOException {
		if (getSettings().hasGeneOntologyAnnotations()) {
			setGoAnnotationReference(ReferenceGoAnnotation.parseGoAnnotationReference(getUniqueBlastResultShortAccessions()));
		}
	}

	public void filterBestScoringBlastResults(Protein prot) {
		for (String blastDatabaseName : prot.getBlastResults().keySet()) {
			prot.getBlastResults().put(blastDatabaseName,
					BlastResult.filterBestScoringBlastResults(prot.getBlastResults().get(blastDatabaseName), 200));
		}
	}

	/**
	 * Method initializes the AHRD-run: 1. Loads Proteins 2. Parses BlastResults
	 * 3. Parses InterproResults 4. Parses database Gene-Ontology-Annotations
	 * 
	 * @throws IOException
	 * @throws MissingAccessionException
	 * @throws MissingProteinException
	 * @throws SAXException
	 * @throws ParsingException
	 */
	public void setup(boolean writeLogMsgs)
			throws IOException, MissingAccessionException, MissingProteinException, SAXException, ParsingException {
		if (writeLogMsgs)
			System.out.println("Started AHRD...\n");

		takeTime();

		initializeProteins();
		if (writeLogMsgs)
			System.out.println("...initialised proteins in " + takeTime() + "sec, currently occupying "
					+ takeMemoryUsage() + " MB");

		// multiple blast-results against different Blast-Databases
		parseBlastResults();
		if (writeLogMsgs)
			System.out.println("...parsed blast results in " + takeTime() + "sec, currently occupying "
					+ takeMemoryUsage() + " MB");

		// GO Annotation Reference (for Proteins in the searched Blast Databases)
		setUpGoAnnotationReference();
		if (writeLogMsgs) {
			System.out.println("...parsed Gene Ontology Annotation (GOA) Reference in " + takeTime()
					+ "sec, currently occupying " + takeMemoryUsage() + " MB");
		}

		// one single InterproResult-File
		if (getSettings().hasValidInterproDatabaseAndResultFile()) {
			InterproResult.initialiseInterproDb();
			parseInterproResult();
			if (writeLogMsgs)
				System.out.println("...parsed interpro results in " + takeTime() + "sec, currently occupying "
						+ takeMemoryUsage() + " MB");
		}
	}

	/**
	 * Assign a HumanReadableDescription to each Protein
	 * 
	 * @throws MissingInterproResultException
	 * @throws IOException
	 * @throws SQLException
	 */
	public void assignHumanReadableDescriptions() throws MissingInterproResultException, IOException, SQLException {
		for (String protAcc : getProteins().keySet()) {
			Protein prot = getProteins().get(protAcc);
			// Find best scoring Blast-Hit's Description-Line (based on
			// evalue):
			filterBestScoringBlastResults(prot);
			// Tokenize each BlastResult's Description-Line and
			// assign the Tokens their Scores:
			// tokenizeBlastResultDescriptionLines(prot);
			prot.getTokenScoreCalculator().assignTokenScores();
			// Tell informative from non-informative Tokens.
			// Assign each non-informative a new Score :=
			// currentScore - (Token-High-Score * Informative-Token-Threshold)
			prot.getTokenScoreCalculator().filterTokenScores();
			// Find the highest scoring Blast-Result:
			prot.getDescriptionScoreCalculator().findHighestScoringBlastResult();
			// filter for each protein's most-informative
			// interpro-results
			InterproResult.filterForMostInforming(prot);
		}
	}
	/**
	 * Assign Gene Ontology terms to each Protein
	 */
	public void assignGeneOntologyTerms() throws MissingInterproResultException, IOException, SQLException {
		// Load a Map of all GO terms
		if (goDB == null) {
			goDB = new GOdatabase().getMap();
		}
		for (Protein protein : this.getProteins().values()) {
			// calculate total and cumulative go term scores 
			Map<String, Double> cumulativeGoTermBitScores = new HashMap<String, Double>();
			Map<String, Double> cumulativeGoTermBlastDatabaseScores = new HashMap<String, Double>();
			Map<String, Double> cumulativeGoTermOverlapScores = new HashMap<String, Double>();
			double totalGoTermBitScore = 0;
			double totalGoTermBlastDatabaseScore = 0;
			double totalGoTermOverlapScore = 0;
			double maxBitScore = 0;
			for (String blastDbName : protein.getBlastResults().keySet()) {
				for (BlastResult blastResult : protein.getBlastResults().get(blastDbName)) {
					totalGoTermBitScore += blastResult.getBitScore();
					totalGoTermBlastDatabaseScore += getSettings().getBlastDbWeight(blastDbName);
					// calculate overlap score
					double overlapScore = TokenScoreCalculator.overlapScore(blastResult.getQueryStart(), blastResult.getQueryEnd(),
							protein.getSequenceLength(), blastResult.getSubjectStart(), blastResult.getSubjectEnd(), blastResult.getSubjectLength());
					totalGoTermOverlapScore += overlapScore; 
					Set<ReferenceGoAnnotation> reference = this.getGoAnnotationReference().get(blastResult.getShortAccession());
					if (reference != null) {
						for (ReferenceGoAnnotation annotation : reference) {
							String goTerm = annotation.getGoTerm();
							// calculate cumulative bit score
							if (!cumulativeGoTermBitScores.containsKey(goTerm)) {
								cumulativeGoTermBitScores.put(goTerm, new Double(blastResult.getBitScore()));
							} else {
								cumulativeGoTermBitScores.put(goTerm, new Double(blastResult.getBitScore() + cumulativeGoTermBitScores.get(goTerm)));
							}
							// calculate cumulative blast database score
							if (!cumulativeGoTermBlastDatabaseScores.containsKey(goTerm)) {
								cumulativeGoTermBlastDatabaseScores.put(goTerm, new Double(getSettings().getBlastDbWeight(blastDbName)));
							} else {
								cumulativeGoTermBlastDatabaseScores.put(goTerm, new Double(getSettings().getBlastDbWeight(blastDbName) + cumulativeGoTermBlastDatabaseScores.get(goTerm)));
							}
							// calculate cumulative overlap score
							if (!cumulativeGoTermOverlapScores.containsKey(goTerm)) {
								cumulativeGoTermOverlapScores.put(goTerm, overlapScore);
							} else {
								cumulativeGoTermOverlapScores.put(goTerm, overlapScore + cumulativeGoTermOverlapScores.get(goTerm));
							}
						}
					}
					// calculate max bit score
					if (blastResult.getBitScore() > maxBitScore) {
						maxBitScore = new Double(blastResult.getBitScore());
					}
				}
			}
			// Calculate GO term evidence code scores
			Map<String, Double> cumulativeGoTermEvidenceCodeWeights = new HashMap<String, Double>();
			Map<String, Integer> termAnnotationCounts = new HashMap<String, Integer>();
			for (String blastDbName : protein.getBlastResults().keySet()) {
				for (BlastResult blastResult : protein.getBlastResults().get(blastDbName)) {
					Set<ReferenceGoAnnotation> reference = this.getGoAnnotationReference().get(blastResult.getShortAccession());
					if (reference != null) {
						for (ReferenceGoAnnotation annotation : reference) {
							String goTerm = annotation.getGoTerm();
							String code = annotation.getEvidenceCode();
							// calculate cumulative evidence code score
							if (!cumulativeGoTermEvidenceCodeWeights.containsKey(goTerm)) {
								cumulativeGoTermEvidenceCodeWeights.put(goTerm, getSettings().getEvidenceCodeWeight(code));
							} else {
								cumulativeGoTermEvidenceCodeWeights.put(goTerm, getSettings().getEvidenceCodeWeight(code) + cumulativeGoTermEvidenceCodeWeights.get(goTerm));
							}
							// calculate number of term annotations 
							if (!termAnnotationCounts.containsKey(goTerm)) {
								termAnnotationCounts.put(goTerm, 1);
							} else {
								termAnnotationCounts.put(goTerm, termAnnotationCounts.get(goTerm) + 1);
							}

						}
					}
				}
			}
			// Calculate GO Term-Scores and GO term high score
			Map<GOterm, Double> goTermScores = new HashMap<GOterm, Double>();
			double goTermHighScore = 0.0;
			for (String blastDbName : protein.getBlastResults().keySet()) {
				for (BlastResult blastResult : protein.getBlastResults().get(blastDbName)) {
					Set<ReferenceGoAnnotation> reference = this.getGoAnnotationReference().get(blastResult.getShortAccession());
					if (reference != null) {
						for (ReferenceGoAnnotation annotation : reference) {
							String termAcc = annotation.getGoTerm();
							GOterm term = goDB.get(termAcc);
							double infoContentScore = 1 - (getSettings().getGoTermScoreInformationContentWeight() * term.getAnnotationFrequency());
							double evidenceCodeScore = 1 - (getSettings().getGoTermScoreEvidenceCodeScoreWeight() * (1 - (cumulativeGoTermEvidenceCodeWeights.get(termAcc) / termAnnotationCounts.get(termAcc)))); 
							double goTermAbundancyScore = getSettings().getTokenScoreBitScoreWeight() * cumulativeGoTermBitScores.get(termAcc) / totalGoTermBitScore 
														+ getSettings().getTokenScoreDatabaseScoreWeight() * cumulativeGoTermBlastDatabaseScores.get(termAcc) / totalGoTermBlastDatabaseScore
														+ getSettings().getTokenScoreOverlapScoreWeight() * cumulativeGoTermOverlapScores.get(termAcc) / totalGoTermOverlapScore;
							double goTermScore = goTermAbundancyScore * infoContentScore * evidenceCodeScore;
							goTermScores.put(term, goTermScore);
							if (goTermScore > goTermHighScore) {
								goTermHighScore = goTermScore;
							}
						}
					}
				}
			}
			// Annotate protein with all goTerms while using goTermScore as annotation confidence
			protein.setGoResultsTermsConfidence(new HashMap<GOterm, Double>(goTermScores));
			// Filter GO Term-Scores
			for (GOterm goTerm : goTermScores.keySet()) {
				if (goTermScores.get(goTerm) < goTermHighScore * getSettings().getInformativeTokenThreshold()) {
					goTermScores.put(goTerm, new Double(goTermScores.get(goTerm) - goTermHighScore * getSettings().getInformativeTokenThreshold()));
				}
			}
			// Find highest scoring GO annotation
			double goAnnotationTopScore = 0.0;
			BlastResult highestScoringBlastResult = null;
			for (String blastDbName : protein.getBlastResults().keySet()) {
				for (BlastResult blastResult : protein.getBlastResults().get(blastDbName)) {
					double sumGoTermScores = 0.0;
					int informativeGoTermCount = 0;
					int goTermCount = 0;
					Set<ReferenceGoAnnotation> reference = this.getGoAnnotationReference().get(blastResult.getShortAccession());
					if (reference != null) { 
						for (ReferenceGoAnnotation annotation : reference) {
							Double goTermScore = goTermScores.get(goDB.get(annotation.getGoTerm()));
							sumGoTermScores += goTermScore * getSettings().getEvidenceCodeWeight(annotation.getEvidenceCode());
							goTermCount++;
							if (goTermScore > goTermHighScore * getSettings().getInformativeTokenThreshold()) {
								informativeGoTermCount++;
							}
						}
					}
					double correctionFactor = ((double) informativeGoTermCount) / ((double) goTermCount);
					double lexicalScore = correctionFactor * sumGoTermScores / goTermHighScore;
					double relativeBlastScore = getSettings().getDescriptionScoreBitScoreWeight(blastDbName) * blastResult.getBitScore() / maxBitScore;
					double goAnnotationScore = lexicalScore + relativeBlastScore;
					if (goAnnotationScore > goAnnotationTopScore) {
						goAnnotationTopScore = goAnnotationScore;
						highestScoringBlastResult = blastResult;
					}
				}
			}
			if (highestScoringBlastResult != null) {
				Set<ReferenceGoAnnotation> reference = getGoAnnotationReference().get(highestScoringBlastResult.getShortAccession());
				Set<String> results = new HashSet<String>();
				for (ReferenceGoAnnotation annotation : reference) {
					results.add(annotation.getGoTerm());
				}
				protein.setGoResults(results);
			} else {
				protein.setGoResults(new HashSet<String>());
			}
		}
	}

	/**
	 * If AHRD is requested to annotate GO terms in accordance to a GO slim set
	 * 
	 * @throws IOException
	 * @throws MissingAccessionException
	 */
	public void annotateWithGoSlim() throws IOException, MissingAccessionException {
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoSlimFile()) {
			Set<GOterm> goSlim = new HashSet<GOterm>();
			// Load a Map of all GO terms
			if (goDB == null) {
				goDB = new GOdatabase().getMap();
			}
			// Load set of GO slim terms
			for (String goSlimFileEntry : getSettings().getGoSlimFile()) {
				Pattern p = Settings.getGoSlimFileGotermRegex();
				Matcher m = p.matcher(goSlimFileEntry);
				if (m.find()) {
					String termAcc = m.group("goTerm");
					GOterm term = goDB.get(termAcc);
					if (term == null)
						throw new MissingAccessionException("Could not find GO term for accession '" + termAcc + "'");
					goSlim.add(term);
				}
			}
			// Annotate proteins:
			// The ancestry of every detailed term is examined for intersections
			// with the GO-Slim set.
			// From the intersection only the term with the highest information
			// content is annotated because some GO-Slim categories exclude some
			// of their child terms.
			for (Iterator<Protein> protIter = getProteins().values().iterator(); protIter.hasNext();) {
				Protein prot = protIter.next();
				for (String termAcc : prot.getGoResults()) {
					GOterm term = goDB.get(termAcc);
					if (term == null) {
						throw new MissingAccessionException("Could not find GO term for accession '" + termAcc + "'");
					}
					Double maxInfoContent = 0.0;
					GOterm highestInfoContentGoSlimTerm = null;
					for (GOterm ancestor : term.getAncestry()) {
						if (goSlim.contains(ancestor) && ancestor.getInformationContent() > maxInfoContent) {
							highestInfoContentGoSlimTerm = ancestor;
						}
					}
					if (highestInfoContentGoSlimTerm != null) {
						prot.getGoSlimTerms().add(highestInfoContentGoSlimTerm);	
					}
				}
			}
		}
	}

	public void termCentricAnnotation() throws IOException, MissingAccessionException {
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoTermCentricTermsFile()) {
			// Load a Map of all GO terms
			if (goDB == null) {
				goDB = new GOdatabase().getMap();
			}
			// Load set of GO centric terms
			for (String goCentricTermFileEntry : getSettings().getGoTermCentricTerms()) {
				Pattern p = Settings.GO_TERM_CENTRIC_TERMS_FILE_GOTERM_REGEX;
				Matcher m = p.matcher(goCentricTermFileEntry);
				if (m.find()) {
					String termAcc = m.group("goTerm");
					GOterm term = goDB.get(termAcc);
					if (term == null)
						throw new MissingAccessionException("Could not find GO term for accession '" + termAcc + "'");
					goCentricTerms.add(term);
				}
			}
			// Determine protein-goTerm association confidence (term centric annotation):
			// The term itself and its child terms are considered for each term.
			// From these terms the maximum confidence is used
			for (GOterm term : this.goCentricTerms) {
				Double maxConfidence = 0.0;
				for (Protein protein : getProteins().values()) {
					Double confidence = 0.0;
					for (GOterm potentialChild : protein.getGoResultsTermsConfidence().keySet()) {
						if (potentialChild.getAncestry().contains(term)) {
							if (protein.getGoResultsTermsConfidence().get(potentialChild) > confidence) {
								confidence = protein.getGoResultsTermsConfidence().get(potentialChild);
							}
						}
					}
					protein.getGoCentricTermConfidences().put(term.getAccession(), confidence);
					// Determine highest association confidence between this GOterm and all proteins  
					if (confidence > maxConfidence) {
						maxConfidence = confidence;
					}
				}
				System.out.println(term.getAccession() + " maxConfidence: " + maxConfidence);
				// Scale the protein-GOterm association confidences according to the highest one so they end up using all the 'space' between 0 and 1 
				for (Protein protein : getProteins().values()) {
					protein.getGoCentricTermConfidences().put(term.getAccession(), protein.getGoCentricTermConfidences().get(term.getAccession()) / maxConfidence);
				}
			}
		}
	}

	public Map<String, Protein> getProteins() {
		return proteins;
	}

	public void setProteins(Map<String, Protein> proteins) {
		this.proteins = proteins;
	}

	public Map<String, Double> getDescriptionScoreBitScoreWeights() {
		return descriptionScoreBitScoreWeights;
	}

	public void setDescriptionScoreBitScoreWeights(Map<String, Double> descriptionScoreBitScoreWeights) {
		this.descriptionScoreBitScoreWeights = descriptionScoreBitScoreWeights;
	}

	public Map<String, Set<ReferenceGoAnnotation>> getGoAnnotationReference() {
		return goAnnotationReference;
	}

	public void setGoAnnotationReference(Map<String, Set<ReferenceGoAnnotation>> goaReference) {
		this.goAnnotationReference = goaReference;
	}

	public Set<String> getUniqueBlastResultShortAccessions() {
		return uniqueBlastResultShortAccessions;
	}

	public void setUniqueBlastResultShortAccessions(Set<String> uniqueBlastResultShortAccessions) {
		this.uniqueBlastResultShortAccessions = uniqueBlastResultShortAccessions;
	}

	public Map<String, GOterm> getGoDB() {
		return goDB;
	}

	public void setGoDB(Map<String, GOterm> goDB) {
		this.goDB = goDB;
	}

}
