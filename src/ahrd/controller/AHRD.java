package ahrd.controller;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Settings.setSettings;
import static ahrd.model.DatabaseGoAnnotations.parseDatabaseGoAnnotations;

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
import ahrd.view.FastaOutputWriter;
import ahrd.view.OutputWriter;
import ahrd.view.TsvOutputWriter;
import nu.xom.ParsingException;

public class AHRD {

	public static final String VERSION = "3.11";

	private Map<String, Protein> proteins;
	private Map<String, Double> descriptionScoreBitScoreWeights = new HashMap<String, Double>();
	private Map<String, Set<String>> databaseGoAnnotations;
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
			// Log
			System.out.println("...assigned highestest scoring human readable descriptions in " + ahrd.takeTime()
					+ "sec, currently occupying " + ahrd.takeMemoryUsage() + " MB");
			// If requested assign GO slim terms in addition to detailed GO term
			// annotation
			ahrd.annotateWithGoSlim();
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
	 * <code>databaseGoAnnotations</code> required only if AHRD is requested to
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
			this.setDatabaseGoAnnotations(new HashMap<String, Set<String>>());
		}
	}

	public void initializeProteins() throws IOException, MissingAccessionException {
		setProteins(Protein.initializeProteins(getSettings().getProteinsFasta()));
	}

	public void parseBlastResults() throws IOException, MissingProteinException, SAXException {
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
	public void setUpDatabaseGoAnnotations() throws IOException {
		if (getSettings().hasGeneOntologyAnnotations()) {
			setDatabaseGoAnnotations(parseDatabaseGoAnnotations(getUniqueBlastResultShortAccessions()));
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

		// Database GO Annotations (for Proteins in the searched Blast Databases)
		setUpDatabaseGoAnnotations();
		if (writeLogMsgs) {
			System.out.println("...parsed database Gene Ontology Annotations (GOA) in " + takeTime()
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
			// currentScore - (Token-High-Score / 2)
			prot.getTokenScoreCalculator().filterTokenScores();
			// Find the highest scoring Blast-Result:
			prot.getDescriptionScoreCalculator().findHighestScoringBlastResult(this.getDatabaseGoAnnotations());
			// If AHRD is requested to annotate Gene Ontology Terms, do so:
			if (getSettings().hasGeneOntologyAnnotations()
					&& prot.getDescriptionScoreCalculator().getHighestScoringBlastResult() != null
					&& getDatabaseGoAnnotations().containsKey(
							prot.getDescriptionScoreCalculator().getHighestScoringBlastResult().getShortAccession())) {
				prot.setGoResults(getDatabaseGoAnnotations()
						.get(prot.getDescriptionScoreCalculator().getHighestScoringBlastResult().getShortAccession()));
			} else {
				prot.setGoResults(new HashSet<String>());
			}
			// filter for each protein's most-informative
			// interpro-results
			InterproResult.filterForMostInforming(prot);
		}
	}

	/**
	 * If AHRD is requested to annotate GO term in accordance to a GO slim set
	 * 
	 * @throws IOException
	 * @throws MissingAccessionException
	 */
	public void annotateWithGoSlim() throws IOException, MissingAccessionException {
		if (getSettings().hasGeneOntologyAnnotations() && getSettings().hasGoSlimFile()) {
			Set<GOterm> goSlim = new HashSet<GOterm>();
			// Load a Map of all GO terms
			goDB = new GOdatabase().getMap();
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

	public Map<String, Set<String>> getDatabaseGoAnnotations() {
		return databaseGoAnnotations;
	}

	public void setDatabaseGoAnnotations(Map<String, Set<String>> databaseGoAnnotations) {
		this.databaseGoAnnotations = databaseGoAnnotations;
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
