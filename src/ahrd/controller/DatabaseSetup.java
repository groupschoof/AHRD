package ahrd.controller;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Settings.setSettings;
import static ahrd.model.AhrdDb.closeDb;
import static ahrd.model.AhrdDb.initializeDb;
import static ahrd.model.ReferenceProtein.parseBlastDatabase;
import static ahrd.model.ReferenceProtein.parseReferenceGoAnnotations;

import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Files;

import com.sleepycat.je.DatabaseException;

public class DatabaseSetup {

	/**
	 * Read input and trigger creation or updating of AHRD's persistent
	 * database.
	 * 
	 * @param args
	 * @throws DatabaseException
	 * @throws IOException
	 */
	public static void main(String[] args) throws DatabaseException, IOException {
		try {
			System.out
					.println("Usage:\njava -Xmx2g -cp ahrd.jar ahrd.controller.DatabaseSetup databaseSetupInput.yml\n");
			setSettings(new Settings(args[0]));
			createOrUpdateAhrdDatabase(true);
			System.out.println("Parsed and persisted reference sequence information as given in '" + args[0]
					+ "'. Database has been saved in '" + getSettings().getAhrd_db() + "'.\n");
		} finally {
			// Close the reference protein persistent store:
			closeDb();
		}
	}

	/**
	 * Checks whether AHRD's database file already exists.
	 * 
	 * @return true if and only if the file as indicated in Settings exists.
	 */
	public static boolean ahrdDatabaseExists() {
		return Files.exists(FileSystems.getDefault().getPath(getSettings().getAhrd_db()));
	}

	/**
	 * Creates or updates AHRD's database. Parses the reference sequence
	 * databases (Fasta-Format) and optionally the Gene Ontology Annotation
	 * (GOA) files, writing parsed information persistently into the Database.
	 * 
	 * @param writeLogMsgs
	 *            - set to true if log messages should be written to standard
	 *            out and error
	 * @throws DatabaseException
	 * @throws IOException
	 */
	public static void createOrUpdateAhrdDatabase(boolean writeLogMsgs) throws DatabaseException, IOException {
		String action = ahrdDatabaseExists() ? "UPDATING EXISTING" : "CREATING NEW";
		if (writeLogMsgs)
			System.out.println(action + " AHRD-Database in directory '"
					+ FileSystems.getDefault().getPath(getSettings().getAhrd_db()).toAbsolutePath() + "'.\n");
		// Open the Database allowing writes ('readOnly' = false)
		initializeDb(false);
		// Parse the reference sequence fasta databases:
		parseReferenceSequenceDatabases(writeLogMsgs);
		// Optionally parse Gene Ontology Annotations (GOA):
		parseGeneOntologyAnnotations(writeLogMsgs);
	}

	/**
	 * Parses the reference sequence databases indicated in the input. Parsed
	 * information is extracted into instances of ReferenceProteins and stored
	 * in AHRD's database persistently.
	 * 
	 * @param writeLogMsgs
	 *            - set to true if log messages should be written to standard
	 *            out and error
	 * @throws IOException
	 */
	public static void parseReferenceSequenceDatabases(boolean writeLogMsgs) throws IOException {
		if (!getSettings().getBlastDatabases().isEmpty()) {
			for (String blastDb : getSettings().getBlastDatabases()) {
				if (writeLogMsgs)
					System.out.println(
							"Starting to parse reference sequence database (Fasta-Format): '" + blastDb + "'.");
				parseBlastDatabase(blastDb);
			}
			if (writeLogMsgs)
				System.out.println("Done parsing reference sequence databases.");
		} else if (writeLogMsgs)
			System.err.println("WARNING: No reference sequence databases (Fasta-Format) were given in input-file.");
	}

	/**
	 * Optionally parses Gene Ontology Annotation (GOA) files. Informs the User
	 * about the process.
	 * 
	 * @param writeLogMsgs
	 *            - set to true if log messages should be written to standard
	 *            out and error
	 * @throws IOException
	 */
	public static void parseGeneOntologyAnnotations(boolean writeLogMsgs) throws IOException {
		if (getSettings().hasGeneOntologyAnnotations()) {
			if (writeLogMsgs)
				System.out.println("Parsing Gene Ontology Annotations (GOA).");
			parseReferenceGoAnnotations();
		} else {
			if (writeLogMsgs)
				System.out.println(
						"No Gene Ontology Annotation (GOA) file was set in input. Skipping parsing of GOA data.");
		}
	}

	/**
	 * Connects to an existing database, if one is found, or creates and feeds a
	 * new one.
	 * 
	 * @throws DatabaseException
	 * @throws IOException
	 */
	public static void setupOrUseExistingDatabase(boolean writeLogMsgs) throws DatabaseException, IOException {
		if (ahrdDatabaseExists()) {
			initializeDb(true);
			if (writeLogMsgs)
				System.out.println("Using EXISTING database in '" + getSettings().getAhrd_db() + "'.");
		} else
			createOrUpdateAhrdDatabase(writeLogMsgs); // User will be informed.
	}
}
