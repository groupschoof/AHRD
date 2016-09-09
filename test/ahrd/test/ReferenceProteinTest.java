package ahrd.test;

import static ahrd.controller.Settings.setSettings;
import static ahrd.model.AhrdDb.closeDb;
import static ahrd.model.AhrdDb.getReferenceProteinDAO;
import static ahrd.model.AhrdDb.initializeDb;
import static ahrd.model.ReferenceProtein.parseBlastDatabase;
import static ahrd.model.ReferenceProtein.parseReferenceGoAnnotations;
import static ahrd.test.TestUtils.initTestSettings;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.Settings;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.ReferenceProtein;
import nu.xom.ParsingException;

public class ReferenceProteinTest {

	/**
	 * Helper method to create ReferenceProteins in AHRD's Database.
	 */
	private void loadReferenceProteins() {
		getReferenceProteinDAO().byAccession.put(new ReferenceProtein("sp|acc_1|accession_one",
				"HumanReadableDescription_1", new Integer(123), "swissprot"));
		getReferenceProteinDAO().byAccession.put(new ReferenceProtein("sp|acc_2|accession_two",
				"HumanReadableDescription_2", new Integer(321), "swissprot"));
	}

	/**
	 * Helper method to check whether the expected ReferenceProteins are in the
	 * (persistent) AHRD-Database.
	 */
	private void testReferenceProteins() {
		ReferenceProtein p1 = getReferenceProteinDAO().byAccession.get("sp|acc_1|accession_one");
		assertNotNull(p1);
		assertEquals("sp|acc_1|accession_one", p1.getAccession());
		assertEquals("HumanReadableDescription_1", p1.getHrd());
		assertEquals(new Integer(123), p1.getSequenceLength());
		assertEquals("swissprot", p1.getBlastDatabaseName());
		assertEquals("acc_1", p1.getShortAccession());
		ReferenceProtein p2 = getReferenceProteinDAO().byAccession.get("sp|acc_2|accession_two");
		assertNotNull(p2);
		assertEquals("sp|acc_2|accession_two", p2.getAccession());
		assertEquals("HumanReadableDescription_2", p2.getHrd());
		assertEquals(new Integer(321), p2.getSequenceLength());
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
			parseBlastDatabase("swissprot");
			assertNotNull(getReferenceProteinDAO().byShortAccession.get("Q9ZWC8"));
			ReferenceProtein rp1 = getReferenceProteinDAO().byAccession.get("sp|Q9LRP3|Y3174_ARATH");
			assertNotNull(rp1);
			assertEquals("Q9LRP3", rp1.getShortAccession());
			assertEquals("Probable receptor-like protein kinase At3g17420", rp1.getHrd());
			assertEquals("swissprot", rp1.getBlastDatabaseName());
			assertEquals(new Integer(467), rp1.getSequenceLength());
			assertNull(getReferenceProteinDAO().byAccession.get("sp|ThisAccessionDoesNotExists|it's_true_believe_me"));
		} finally {
			closeDb();
		}
	}

	@Test
	public void testParseReferenceGoAnnotations()
			throws IOException, MissingAccessionException, MissingProteinException, SAXException, ParsingException {
		setSettings(new Settings("./test/resources/ahrd_input_seq_sim_table_go_prediction.yml"));
		try {
			initializeDb(false);
			parseBlastDatabase("swissprot");
			parseBlastDatabase("trembl");
			parseBlastDatabase("tair");
			parseReferenceGoAnnotations();
			ReferenceProtein rp1 = getReferenceProteinDAO().byShortAccession.get("W9QFR0");
			assertNotNull("ReferenceProtein 'W9QFR0' should be in the Database, but is null!", rp1);
			assertNotNull("ReferenceProtein 'W9QFR0' should have GO terms assigned, but its GO term set is null!",
					rp1.getGoTerms());
			assertTrue("ReferenceProtein 'W9QFR0' should have GO terms assigned, but its GO term set is empty!",
					!rp1.getGoTerms().isEmpty());
			assertTrue("ReferenceProtein 'W9QFR0' should have GO term 'GO:0009058' assigned, but has not!",
					rp1.getGoTerms().contains("GO:0009058"));
			assertTrue("ReferenceProtein 'W9QFR0' should have GO term 'GO:0030170' assigned, but has not!",
					rp1.getGoTerms().contains("GO:0030170"));
			ReferenceProtein rp2 = getReferenceProteinDAO().byShortAccession.get("AT1G01040.1");
			assertNotNull("ReferenceProtein 'AT1G01040.1' should be in the Database, but is null!", rp2);
			assertNotNull("ReferenceProtein 'AT1G01040.1' should have GO terms assigned, but its GO term set is null!",
					rp2.getGoTerms());
			assertTrue("ReferenceProtein 'AT1G01040.1' should have GO terms assigned, but its GO term set is empty!",
					!rp2.getGoTerms().isEmpty());
			assertTrue("ReferenceProtein 'AT1G01040.1' should have GO term 'GO:0003824' assigned, but has not!",
					rp2.getGoTerms().contains("GO:0003824"));
			assertTrue("ReferenceProtein 'AT1G01040.1' should have GO term 'GO:0003870' assigned, but has not!",
					rp2.getGoTerms().contains("GO:0003870"));
		} finally {
			closeDb();
		}
	}
}
