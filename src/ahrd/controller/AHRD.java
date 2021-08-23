package ahrd.controller;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Settings.setSettings;
import static ahrd.controller.Settings.SHORT_ACCESSION_GROUP_NAME;
import static ahrd.controller.Settings.DEFAULT_SHORT_ACCESSION_REGEX;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.logging.Handler;
import java.util.logging.Logger;

import org.semanticweb.owlapi.model.OWLOntologyCreationException;

import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.BlastResult;
import ahrd.model.GOdatabase;
import ahrd.model.GOterm;
import ahrd.model.Protein;
import ahrd.model.TokenScoreCalculator;
import ahrd.model.ReferenceGoAnnotation;
import ahrd.view.FastaOutputWriter;
import ahrd.view.OutputWriter;
import ahrd.view.TsvOutputWriter;
import ahrd.view.DualConsoleHandler;

public class AHRD {

	public static final String VERSION = "3.11";
	
	protected static final Logger LOGGER = Logger.getLogger("global");
	
	private Map<String, Protein> proteins;
	private Map<String, Protein> shortAccsProteins;
	private Map<String, Set<ReferenceGoAnnotation>> goAnnotationReference;
	private Set<String> uniqueBlastResultShortAccessions;
	private long timestamp;
	private long memorystamp;
	protected Map<String, GOterm> goDB;

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
		System.out.println("Usage:\njava -jar ahrd.jar input.yml\n");

		try {
			AHRD ahrd = new AHRD(args[0]);
			// Load and parse all inputs
			ahrd.setup();
			// After the setup the unique short accessions are no longer needed:
			ahrd.setUniqueBlastResultShortAccessions(null);

			// Iterate over all Proteins and assign the best scoring Human
			// Readable Description
			ahrd.assignHumanReadableDescriptions();
			// Log
			LOGGER.info("Assigned highestest scoring human readable descriptions in " + ahrd.takeTime() + "sec, currently occupying " + ahrd.takeMemoryUsage() + " MB");
			// If requested iterate over all Proteins and assign the best scoring Gene Ontology terms
			if (getSettings().doAnnotateGoTerms()) {
				ahrd.assignGeneOntologyTerms();
				LOGGER.info("Assigned highestest scoring GO terms in " + ahrd.takeTime() + "sec, currently occupying " + ahrd.takeMemoryUsage() + " MB");
			}
			// If requested assign GO slim terms in addition to detailed GO term
			// annotation
			ahrd.annotateWithGoSlim();
			// Write result to output-file:
			LOGGER.info("Writing output to '" + getSettings().getPathToOutput() + "'.");
			OutputWriter ow = initializeOutputWriter(ahrd.getProteins().values());
			ow.writeOutput();
			LOGGER.info("Wrote output in " + ahrd.takeTime() + "sec, currently occupying " + ahrd.takeMemoryUsage() + " MB");
		} catch (Exception e) {
			LOGGER.severe("We are sorry, an un-expected ERROR occurred:");
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
		// Replace all previous logging handlers with one that pushes SEVERE and WARNING to stderr and all other levels to stdout 
		LOGGER.setUseParentHandlers(false);
		for(Handler handler : LOGGER.getHandlers()) {
		    LOGGER.removeHandler(handler);
		}
		LOGGER.addHandler(new DualConsoleHandler());
		// Parse the settings from the YML input
		setSettings(new Settings(pathToYmlInput));
		// Limit the size of Java's default common threadpool to the specified number of threads
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(getSettings().getNthreads() - 1));
		// If GO annotation references are provided for any blast database GO annotation can be performed
		for (String blastDatabaseName : getSettings().getBlastDatabases()) {
			if (getSettings().getPathToGeneOntologyReference(blastDatabaseName) != null
					&& new File(getSettings().getPathToGeneOntologyReference(blastDatabaseName)).exists()) {
				getSettings().setAnnotateGoTerms(true);
			}
		}
		// The following fields are only used if AHRD is requested to generate
		// Gene Ontology term annotations:
		if (getSettings().doAnnotateGoTerms()) {
			this.setUniqueBlastResultShortAccessions(new HashSet<String>());
			this.setGoAnnotationReference(new HashMap<String, Set<ReferenceGoAnnotation>>());
		}
	}

	public void initializeProteins() throws IOException, MissingAccessionException {
		setProteins(Protein.initializeProteins(getSettings().getProteinsFasta()));
	}

	public void parseBlastResults() throws IOException, MissingProteinException, MissingAccessionException {
		for (String blastDatabase : getSettings().getBlastDatabases()) {
			BlastResult.readBlastResults(getProteins(), blastDatabase, getUniqueBlastResultShortAccessions());
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
	 * 3. Parses database Gene-Ontology-Annotations
	 * 
	 * @throws IOException
	 * @throws MissingAccessionException
	 * @throws MissingProteinException
	 */
	public void setup()	throws IOException, MissingAccessionException, MissingProteinException {
		LOGGER.info("Started AHRD...");
		takeTime();
		
		initializeProteins();
		LOGGER.info("Initialised proteins in " + takeTime() + "sec, currently occupying " + takeMemoryUsage() + " MB");
		
		// multiple blast-results against different Blast-Databases
		parseBlastResults();
		LOGGER.info("Parsed blast results in " + takeTime() + "sec, currently occupying " + takeMemoryUsage() + " MB");

		// GO Annotation Reference (for Proteins in the searched Blast Databases)
		if (getSettings().doAnnotateGoTerms()) {
			setGoAnnotationReference(ReferenceGoAnnotation.parseGoAnnotationReference(getUniqueBlastResultShortAccessions()));
			LOGGER.info("Parsed Gene Ontology Annotation (GOA) Reference in " + takeTime() + "sec, currently occupying " + takeMemoryUsage() + " MB");
		}
	}

	/**
	 * Assign a HumanReadableDescription to each Protein
	 * 
	 */
	public void assignHumanReadableDescriptions() {
		getProteins().values().parallelStream().forEach(prot -> {
			// Find best scoring Blast-Hit's Description-Line (based on evalue):
			filterBestScoringBlastResults(prot);
			// Tokenize each BlastResult's Description-Line and assign the Tokens their Scores:
			prot.getTokenScoreCalculator().assignTokenScores();
			// Tell informative from non-informative Tokens.
			// Assign each non-informative a new Score :=
			// currentScore - (Token-High-Score * Informative-Token-Threshold)
			prot.getTokenScoreCalculator().filterTokenScores();
			// Find the highest scoring Blast-Result:
			prot.getDescriptionScoreCalculator().findHighestScoringBlastResult();
		});
	}
	
	/**
	 * Assign Gene Ontology terms to each Protein
	 */
	public void assignGeneOntologyTerms() {
		getProteins().values().parallelStream().forEach(protein -> {
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
					totalGoTermBlastDatabaseScore += getSettings().getGoBlastDbWeight(blastDbName);
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
								cumulativeGoTermBitScores.put(goTerm, blastResult.getBitScore());
							} else {
								cumulativeGoTermBitScores.put(goTerm, blastResult.getBitScore() + cumulativeGoTermBitScores.get(goTerm));
							}
							// calculate cumulative blast database score
							if (!cumulativeGoTermBlastDatabaseScores.containsKey(goTerm)) {
								cumulativeGoTermBlastDatabaseScores.put(goTerm, Double.valueOf(getSettings().getGoBlastDbWeight(blastDbName)));
							} else {
								cumulativeGoTermBlastDatabaseScores.put(goTerm, Double.valueOf(getSettings().getGoBlastDbWeight(blastDbName) + cumulativeGoTermBlastDatabaseScores.get(goTerm)));
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
						maxBitScore = blastResult.getBitScore();
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
			Map<String, Double> goTermScores = new HashMap<String, Double>();
			double goTermHighScore = 0.0;
			for (String blastDbName : protein.getBlastResults().keySet()) {
				for (BlastResult blastResult : protein.getBlastResults().get(blastDbName)) {
					Set<ReferenceGoAnnotation> reference = this.getGoAnnotationReference().get(blastResult.getShortAccession());
					if (reference != null) {
						for (ReferenceGoAnnotation annotation : reference) {
							String termAcc = annotation.getGoTerm();
							double evidenceCodeScore = 1 - (getSettings().getGoTermScoreEvidenceCodeScoreWeight() * (1 - (cumulativeGoTermEvidenceCodeWeights.get(termAcc) / termAnnotationCounts.get(termAcc)))); 
							double goTermAbundancyScore = getSettings().getGoTokenScoreBitScoreWeight() * cumulativeGoTermBitScores.get(termAcc) / totalGoTermBitScore 
														+ getSettings().getGoTokenScoreDatabaseScoreWeight() * cumulativeGoTermBlastDatabaseScores.get(termAcc) / totalGoTermBlastDatabaseScore
														+ getSettings().getGoTokenScoreOverlapScoreWeight() * cumulativeGoTermOverlapScores.get(termAcc) / totalGoTermOverlapScore;
							double goTermScore = goTermAbundancyScore * evidenceCodeScore;
							goTermScores.put(termAcc, goTermScore);
							if (goTermScore > goTermHighScore) {
								goTermHighScore = goTermScore;
							}
						}
					}
				}
			}
			// Filter GO Term-Scores
			for (String goTerm : goTermScores.keySet()) {
				if (goTermScores.get(goTerm) < goTermHighScore * getSettings().getGoInformativeTokenThreshold()) {
					goTermScores.put(goTerm, goTermScores.get(goTerm) - goTermHighScore * getSettings().getGoInformativeTokenThreshold());
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
							Double goTermScore = goTermScores.get(annotation.getGoTerm());
							sumGoTermScores += goTermScore * getSettings().getEvidenceCodeWeight(annotation.getEvidenceCode());
							goTermCount++;
							if (goTermScore > goTermHighScore * getSettings().getGoInformativeTokenThreshold()) {
								informativeGoTermCount++;
							}
						}
					}
					double correctionFactor = ((double) informativeGoTermCount) / ((double) goTermCount);
					double lexicalScore = correctionFactor * sumGoTermScores / goTermHighScore;
					double relativeBlastScore = getSettings().getGoScoreBitScoreWeight(blastDbName) * blastResult.getBitScore() / maxBitScore;
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
		});
	}
	
	/**
	 * If AHRD is requested to annotate GO terms in accordance to a GO slim set
	 * @throws IOException 
	 * @throws OWLOntologyCreationException 
	 * @throws MissingAccessionException 
	 * 
	 */
	public void annotateWithGoSlim() throws OWLOntologyCreationException, IOException, MissingAccessionException {
		if (getSettings().doAnnotateGoTerms() && getSettings().hasGoSlimFile()) {
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
							maxInfoContent = ancestor.getInformationContent();
						}
					}
					if (highestInfoContentGoSlimTerm != null) {
						prot.getGoSlimTerms().add(highestInfoContentGoSlimTerm);	
					}
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

	public Map<String, Protein> getShortAccsProteins() {
		if (shortAccsProteins == null) {
			shortAccsProteins = new HashMap<>();
			for (Protein p : getProteins().values()) {
				Matcher m = DEFAULT_SHORT_ACCESSION_REGEX.matcher(p.getAccession());
				String shortAccession = p.getAccession();
				if (!m.find()) {
					LOGGER.warning("Regular Expression '" + DEFAULT_SHORT_ACCESSION_REGEX.toString()
					+ "' does NOT match - using pattern.find(...) - Protein Accession '" + p.getAccession()
					+ "' - continuing with the original accession. This might lead to unrecognized ground thruth GO annotations!");
				} else {
					shortAccession = m.group(SHORT_ACCESSION_GROUP_NAME);
				}
				shortAccsProteins.put(shortAccession, p);
			}
		}
		return shortAccsProteins;
	}

	public void setShortAccsProteins(Map<String, Protein> shortAccsProteins) {
		this.shortAccsProteins = shortAccsProteins;
	}

}
