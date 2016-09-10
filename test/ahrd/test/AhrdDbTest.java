package ahrd.test;

import static ahrd.model.AhrdDb.deleteDb;
import static ahrd.model.AhrdDb.initializeDb;
import static ahrd.model.AhrdDb.getAhrdDbEnv;

import static ahrd.test.TestUtils.initTestSettings;
import static org.junit.Assert.assertEquals;
import static ahrd.controller.Settings.getSettings;

import java.io.IOException;

import org.junit.Test;

import com.sleepycat.je.DatabaseException;
import com.sleepycat.je.EnvironmentMutableConfig;

public class AhrdDbTest {

	@Test
	public void testEnvironment() throws DatabaseException, IOException {
		try {
			initTestSettings();
			initializeDb(false);
			EnvironmentMutableConfig mutableConfig = getAhrdDbEnv().getMutableConfig();
			assertEquals(getSettings().getAhrdDbCachePercent(), mutableConfig.getCachePercent());
		} finally {
			deleteDb();
		}
	}
}
