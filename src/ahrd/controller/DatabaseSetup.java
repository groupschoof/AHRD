package ahrd.controller;

import static ahrd.model.AhrdDb.close;
import static ahrd.model.AhrdDb.initialize;

public class DatabaseSetup {

	public static void main(String[] args) {
		System.out.println("Usage:\njava -Xmx2g -cp ahrd.jar ahrd.controller.DatabaseSetup databaseSetupInput.yml\n");
		try {
			initialize(false);
		} finally {
			// Close the reference protein persistent store:
			close();
		}
	}

}
