package ahrd.controller;

import static ahrd.controller.DatabaseSetup.setupOrUseExistingDatabase;
import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Settings.setSettings;
import static ahrd.model.AhrdDb.closeDb;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.xml.sax.SAXException;

import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingInterproResultException;
import ahrd.exception.MissingProteinException;
import ahrd.model.BlastResult;
import ahrd.model.InterproResult;
import ahrd.model.Protein;
import ahrd.view.FastaOutputWriter;
import ahrd.view.IOutputWriter;
import ahrd.view.OutputWriter;
import nu.xom.ParsingException;

public class AHRD {

	public static final String VERSION = "3.4";

	private Map<String, Protein> proteins;
	private Map<String, Double> descriptionScoreBitScoreWeights = new HashMap<String, Double>();
	private long timestamp;
	private long memorystamp;

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
		System.out.println("Usage:\njava -Xmx30g -jar ahrd.jar input.yml\n");

		try {
			AHRD ahrd = new AHRD(args[0]);
			// Load and parse all inputs
			ahrd.setup(true);

			// Iterate over all Proteins and assign the best scoring Human
			// Readable Description
			ahrd.assignHumanReadableDescriptions();
			// Log
			System.out.println("...assigned highestest scoring human readable descriptions in " + ahrd.takeTime()
					+ "sec, currently occupying " + ahrd.takeMemoryUsage() + " MB");
			// Write result to output-file:
			System.out.println("Writing output to '" + getSettings().getPathToOutput() + "'.");
			IOutputWriter ow = initializeOutputWriter(ahrd.getProteins().values());
			ow.writeOutput();
			// Log
			System.out.println("Wrote output in " + ahrd.takeTime() + "sec, currently occupying "
					+ ahrd.takeMemoryUsage() + " MB");

			System.out.println("\n\nDONE");
		} catch (Exception e) {
			System.err.println("We are sorry, an un-expected ERROR occurred:");
			e.printStackTrace(System.err);
		} finally {
			closeDb();
		}
	}

	public static IOutputWriter initializeOutputWriter(Collection<Protein> proteins) {
		IOutputWriter ow = null;
		if (getSettings().doOutputFasta())
			ow = new FastaOutputWriter(proteins);
		else
			ow = new OutputWriter(proteins);
		return ow;
	}

	/**
	 * Constructor initializes this run's settings as a thread-local variable.
	 * Also conditionally initializes fields
	 * <code>uniqueBlastResultShortAccessions</code> and
	 * <code>referenceGoAnnotations</code> required only if AHRD is requested to
	 * generate Gene Ontology term annotations.
	 * 
	 * @param pathToYmlInput
	 * @throws IOException
	 */
	public AHRD(String pathToYmlInput) throws IOException {
		super();
		setSettings(new Settings(pathToYmlInput));
	}

	public void initializeProteins() throws IOException, MissingAccessionException {
		setProteins(Protein.initializeProteins(getSettings().getProteinsFasta()));
	}

	public void parseBlastResults()
			throws IOException, MissingProteinException, SAXException, MissingAccessionException {
		for (String blastDatabase : getSettings().getBlastDatabases()) {
			BlastResult.readBlastResults(getProteins(), blastDatabase);
		}
	}

	@Deprecated
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

	public void filterBestScoringBlastResults(Protein prot) {
		for (String blastDatabaseName : prot.getBlastResults().keySet()) {
			prot.getBlastResults().put(blastDatabaseName,
					BlastResult.filterBestScoringBlastResults(prot.getBlastResults().get(blastDatabaseName), 200));
		}
	}

	/**
	 * Method initializes the AHRD-run: 1. Loads Proteins 2. Parses BlastResults
	 * 3. Parses InterproResults. Step 3 is deprecated!
	 * 
	 * @throws IOException
	 * @throws MissingAccessionException
	 * @throws MissingProteinException
	 * @throws SAXException
	 * @throws ParsingException
	 */
	public void setup(boolean writeLogMsgs)
			throws IOException, MissingAccessionException, MissingProteinException, SAXException, ParsingException {
		setupOrUseExistingDatabase(writeLogMsgs);

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

		// one single InterproResult-File (DEPRECATED)
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
			prot.getDescriptionScoreCalculator().findHighestScoringBlastResult();
			// filter for each protein's most-informative
			// interpro-results (DEPRECATED)
			InterproResult.filterForMostInforming(prot);
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
}
