package ahrd.test;

import static ahrd.model.AhrdDb.close;
import static ahrd.model.AhrdDb.getReferenceProteinDAO;
import static ahrd.model.AhrdDb.initialize;
import static ahrd.test.TestUtils.initTestSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import org.junit.Test;

import ahrd.model.ReferenceProtein;

public class ReferenceProteinTest {

	private void loadReferenceProteins() {
		getReferenceProteinDAO().byAccession.put(new ReferenceProtein("sp|acc_1|accession_one",
				"HumanReadableDescription_1", new Long(123), "swissprot"));
		getReferenceProteinDAO().byAccession.put(new ReferenceProtein("sp|acc_2|accession_two",
				"HumanReadableDescription_2", new Long(321), "swissprot"));
	}

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

	@Test
	public void testLoadReferenceProteins() throws IOException {
		initTestSettings();
		try {
			initialize(false);
			// Create ReferenceProteins in AHRD's DB:
			loadReferenceProteins();
			// Test their existence:
			testReferenceProteins();
			// Test their PERSISTENT existence:
			close();
			initialize(false);
			testReferenceProteins();
		} finally {
			close();
		}
	}

}
