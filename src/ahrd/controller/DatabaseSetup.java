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
		System.out.println("Usage:\njava -Xmx2g -cp ahrd.jar ahrd.controller.DatabaseSetup databaseSetupInput.yml\n");
		setSettings(new Settings(args[0]));
		craeteOrUpdateAhrdDatabase();
		System.out.println("Parsed and persisted reference sequence information as given in '" + args[0]
				+ "'. Database has been saved in '" + getSettings().getAhrd_db() + "'.\n");
	}

	/**
	 * Creates or updates AHRD's database. Parses the reference sequence
	 * databases (Fasta-Format) and optionally the Gene Ontology Annotation
	 * (GOA) files, writing parsed information persistently into the Database.
	 * 
	 * @throws DatabaseException
	 * @throws IOException
	 */
	public static void craeteOrUpdateAhrdDatabase() throws DatabaseException, IOException {
		String action = Files.exists(FileSystems.getDefault().getPath(getSettings().getAhrd_db())) ? "UPDATING EXISTING"
				: "CREATING NEW";
		System.out.println(action + " AHRD-Database in directory '"
				+ FileSystems.getDefault().getPath(getSettings().getAhrd_db()).toAbsolutePath() + "'.\n");
		try {
			// Open the Database allowing writes ('readOnly' = false)
			initializeDb(false);
			// Parse the reference sequence fasta databases:
			parseReferenceSequenceDatabases();
			// Optionally parse Gene Ontology Annotations (GOA):
			parseGeneOntologyAnnotations();
		} finally {
			// Close the reference protein persistent store:
			closeDb();
		}
	}

	/**
	 * Parses the reference sequence databases indicated in the input. Parsed
	 * information is extracted into instances of ReferenceProteins and stored
	 * in AHRD's database persistently.
	 * 
	 * @throws IOException
	 */
	public static void parseReferenceSequenceDatabases() throws IOException {
		if (!getSettings().getBlastDatabases().isEmpty()) {
			for (String blastDb : getSettings().getBlastDatabases()) {
				System.out.println("Starting to parse reference sequence database (Fasta-Format): '" + blastDb + "'.");
				parseBlastDatabase(blastDb);
			}
			System.out.println("Done parsing reference sequence databases.");
		} else
			System.err.println("WARNING: No reference sequence databases (Fasta-Format) were given in input-file.");
	}

	/**
	 * Optionally parses Gene Ontology Annotation (GOA) files. Informs the User
	 * about the process.
	 * 
	 * @throws IOException
	 */
	public static void parseGeneOntologyAnnotations() throws IOException {
		if (getSettings().hasGeneOntologyAnnotations()) {
			System.out.println("Parsing Gene Ontology Annotations (GOA).");
			parseReferenceGoAnnotations();
		} else {
			System.out
					.println("No Gene Ontology Annotation (GOA) file was set in input. Skipping parsing of GOA data.");
		}
	}
}
