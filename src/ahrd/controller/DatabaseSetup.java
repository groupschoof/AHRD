package ahrd.controller;

import static ahrd.model.AhrdDb.closeDb;
import static ahrd.model.AhrdDb.initializeDb;

import java.io.IOException;

import com.sleepycat.je.DatabaseException;

public class DatabaseSetup {

	public static void main(String[] args) throws DatabaseException, IOException {
		System.out.println("Usage:\njava -Xmx2g -cp ahrd.jar ahrd.controller.DatabaseSetup databaseSetupInput.yml\n");
		try {
			initializeDb(false);
		} finally {
			// Close the reference protein persistent store:
			closeDb();
		}
	}

}
