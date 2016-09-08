package ahrd.test;

import static ahrd.model.AhrdDb.closeDb;
import static ahrd.model.AhrdDb.getReferenceProteinDAO;
import static ahrd.model.AhrdDb.initializeDb;
import static ahrd.test.TestUtils.initTestSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import org.junit.Test;

import ahrd.exception.MissingProteinException;
import ahrd.model.ReferenceProtein;

public class ReferenceProteinTest {

	/**
	 * Helper method to create ReferenceProteins in AHRD's Database.
	 */
	private void loadReferenceProteins() {
		getReferenceProteinDAO().byAccession.put(new ReferenceProtein("sp|acc_1|accession_one",
				"HumanReadableDescription_1", new Long(123), "swissprot"));
		getReferenceProteinDAO().byAccession.put(new ReferenceProtein("sp|acc_2|accession_two",
				"HumanReadableDescription_2", new Long(321), "swissprot"));
	}

	/**
	 * Helper method to check wether the expected ReferenceProteins are in the
	 * (persistent) AHRD-Database.
	 */
	private void testReferenceProteins() {
		ReferenceProtein p1 = getReferenceProteinDAO().byAccession.get("sp|acc_1|accession_one");
		assertNotNull(p1);
		assertEquals("sp|acc_1|accession_one", p1.getAccession());
		assertEquals("HumanReadableDescription_1", p1.getHrd());
		assertEquals(new Long(123), p1.getSequenceLength());
		assertEquals("swissprot", p1.getBlastDatabaseName());
		assertEquals("acc_1", p1.getShortAccession());
		ReferenceProtein p2 = getReferenceProteinDAO().byAccession.get("sp|acc_2|accession_two");
		assertNotNull(p2);
		assertEquals("sp|acc_2|accession_two", p2.getAccession());
		assertEquals("HumanReadableDescription_2", p2.getHrd());
		assertEquals(new Long(321), p2.getSequenceLength());
		assertEquals("swissprot", p2.getBlastDatabaseName());
		assertEquals("acc_2", p2.getShortAccession());
		assertNull(getReferenceProteinDAO().byAccession.get("sp|ThisAccessionDoesNotExists|it's_true_believe_me"));
		assertEquals(p1.getShortAccession(),
				getReferenceProteinDAO().byShortAccession.get("acc_1").getShortAccession());
	}

	/**
	 * Test creation, retrieval and persistence of ReferenceProeins with AHRD's
	 * Database.
	 * 
	 * @throws IOException
	 */
	@Test
	public void testLoadReferenceProteins() throws IOException {
		initTestSettings();
		try {
			initializeDb(false);
			// Create ReferenceProteins in AHRD's DB:
			loadReferenceProteins();
			// Test their existence:
			testReferenceProteins();
			// Test their PERSISTENT existence:
			closeDb();
			initializeDb(true);
			testReferenceProteins();
		} finally {
			closeDb();
		}
	}

	@Test
	public void testParseBlastDatabase() throws IOException, MissingProteinException {
		initTestSettings();
		try {
			initializeDb(false);
			ReferenceProtein.parseBlastDatabase("swissprot");
			assertNotNull(getReferenceProteinDAO().byShortAccession.get("Q9ZWC8"));
			ReferenceProtein rp1 = getReferenceProteinDAO().byAccession.get("sp|Q9LRP3|Y3174_ARATH");
			assertNotNull(rp1);
			assertEquals("Q9LRP3", rp1.getShortAccession());
			assertEquals("Probable receptor-like protein kinase At3g17420", rp1.getHrd());
			assertEquals("swissprot", rp1.getBlastDatabaseName());
			assertEquals(new Long(467), rp1.getSequenceLength());
			assertNull(getReferenceProteinDAO().byAccession.get("sp|ThisAccessionDoesNotExists|it's_true_believe_me"));
		} finally {
			closeDb();
		}
	}

}
