package ahrd.model;

import static ahrd.controller.Settings.getSettings;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;

import com.sleepycat.je.DatabaseException;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.persist.EntityStore;
import com.sleepycat.persist.PrimaryIndex;
import com.sleepycat.persist.SecondaryIndex;
import com.sleepycat.persist.StoreConfig;

public class AhrdDb {

	/**
	 * An accessor to handle storage and retrieval of ReferenceProteins with
	 * Berkeley-DB.
	 */
	public static class ReferenceProteinAccessor {
		public PrimaryIndex<String, ReferenceProtein> byAccession;
		public SecondaryIndex<String, String, ReferenceProtein> byShortAccession;

		public ReferenceProteinAccessor(EntityStore store) {
			byAccession = store.getPrimaryIndex(String.class, ReferenceProtein.class);
			byShortAccession = store.getSecondaryIndex(byAccession, String.class, "shortAccession");
		}
	}

	/*
	 * Thread-Local Variables handling persitent storage
	 */
	private static final ThreadLocal<Environment> ahrdDbEnv = new ThreadLocal<Environment>();
	private static final ThreadLocal<EntityStore> ahrdStore = new ThreadLocal<EntityStore>();
	private static final ThreadLocal<ReferenceProteinAccessor> referenceProteinDAO = new ThreadLocal<ReferenceProteinAccessor>();

	public static ReferenceProteinAccessor getReferenceProteinDAO() {
		return referenceProteinDAO.get();
	}

	public static Environment getAhrdDbEnv() {
		return ahrdDbEnv.get();
	}
	/*
	 * END Thread-Local-Variables
	 */

	/**
	 * Initializes and creates, if necessary, the Database. Stores the Database
	 * and the Environment in a thread-local variable.
	 * 
	 * @param readonly
	 *            Set to true if and only if pure read-access is wanted. This
	 *            should be the case if AHRD is run after Database-Setup has
	 *            been executed.
	 * @throws DatabaseException
	 * @throws IOException
	 */
	public static void initializeDb(boolean readonly) throws DatabaseException, IOException {
		EnvironmentConfig envConfig = new EnvironmentConfig();
		envConfig.setTransactional(false);
		envConfig.setAllowCreate(!readonly);
		envConfig.setReadOnly(readonly);
		envConfig.setCachePercent(getSettings().getAhrdDbCachePercent());
		envConfig.setConfigParam("je.log.fileMax", "1073741824");
		File ahrdDbFile = new File(getSettings().getAhrd_db());
		if (!ahrdDbFile.exists())
			ahrdDbFile.mkdirs();
		ahrdDbEnv.set(new Environment(ahrdDbFile, envConfig));
		StoreConfig storeConfig = new StoreConfig();
		storeConfig.setTransactional(false);
		storeConfig.setAllowCreate(!readonly);
		storeConfig.setReadOnly(readonly);
		storeConfig.setDeferredWrite(!readonly);
		ahrdStore.set(new EntityStore(ahrdDbEnv.get(), "AhrdStore", storeConfig));
		referenceProteinDAO.set(new ReferenceProteinAccessor(ahrdStore.get()));
	}

	/**
	 * Cleans up and closes the connection to the Database-File.
	 */
	public static void closeDb() {
		if (ahrdStore.get() != null) {
			ahrdStore.get().close();
			ahrdStore.set(null);
		}
		if (ahrdDbEnv.get() != null) {
			ahrdDbEnv.get().close();
			ahrdDbEnv.set(null);
		}
		if (referenceProteinDAO.get() != null) {
			referenceProteinDAO.set(null);
		}
	}

	/**
	 * Closes all connections to the database and deletes it from the file
	 * system.
	 * 
	 * @throws IOException
	 */
	public static void deleteDb() throws IOException {
		closeDb();
		Files.walkFileTree(FileSystems.getDefault().getPath(getSettings().getAhrd_db()), new SimpleFileVisitor<Path>() {
			@Override
			public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
				Files.delete(file);
				return FileVisitResult.CONTINUE;
			}

			@Override
			public FileVisitResult postVisitDirectory(Path dir, IOException e) throws IOException {
				if (e == null) {
					Files.delete(dir);
					return FileVisitResult.CONTINUE;
				} else {
					// directory iteration failed
					throw e;
				}
			}
		});
	}
}
