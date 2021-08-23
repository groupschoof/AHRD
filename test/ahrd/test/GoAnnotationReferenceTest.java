package ahrd.test;

import static ahrd.controller.Settings.getSettings;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;

import org.junit.Before;
import org.junit.Test;

import ahrd.controller.AHRD;
import ahrd.controller.Utils;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;
import ahrd.model.Protein;
import ahrd.model.ReferenceGoAnnotation;

public class GoAnnotationReferenceTest {

	private AHRD ahrd;

	@Before
	public void setUp() throws IOException {
		ahrd = new AHRD(
				"./test/resources/ahrd_input_seq_sim_table_go_prediction.yml");
		getSettings().setLoggingLevel(Level.WARNING);
	}

	@Test
	public void testdoAnnotateGoTerms() {
		// Should also predict GO annotations with default test-Settings:
		assertTrue(getSettings().doAnnotateGoTerms());
	}

	@Test
	public void testUniqueShortAccessions() throws IOException, MissingAccessionException, MissingProteinException {
		ahrd.setup();
		assertNotNull(ahrd.getUniqueBlastResultShortAccessions());
		// Somehow assertEquals does not work on Collections as expected, hence
		// the following work-around:
		List<String> expectedUniqueShortAccessions = Utils
				.fromFile("./test/resources/all_blast_hits_blast8_tabular_searches.txt");
		for (String sa : expectedUniqueShortAccessions) {
			assertTrue(ahrd.getUniqueBlastResultShortAccessions().contains(sa));
		}
		assertEquals(expectedUniqueShortAccessions.size(), ahrd
				.getUniqueBlastResultShortAccessions().size());
		ahrd.getUniqueBlastResultShortAccessions().removeAll(
				expectedUniqueShortAccessions);
		assertEquals(0, ahrd.getUniqueBlastResultShortAccessions().size());
	}

	@Test
	public void testParseGoAnnotationReference() throws IOException, MissingAccessionException, MissingProteinException {
		ahrd.setup();
		assertNotNull(ahrd.getGoAnnotationReference());
		assertTrue(!ahrd.getGoAnnotationReference().isEmpty());
		assertEquals(4, ahrd.getGoAnnotationReference().size());
		Set<ReferenceGoAnnotation> reference = ahrd.getGoAnnotationReference().get("AT1G01040");
		Set<String> dbGos = new HashSet<String>();
		for (ReferenceGoAnnotation annotation : reference) {
			dbGos.add(annotation.getGoTerm());
		}
		assertTrue(dbGos.contains("GO:0005634"));
		assertTrue(dbGos.contains("GO:0008026"));
	}

	@Test
	public void testAnnotatesGoTerms() throws IOException, MissingAccessionException, MissingProteinException {
		ahrd.setup();
		ahrd.assignGeneOntologyTerms();
		Protein p = ahrd.getProteins().get("gene:chr01.1056:mRNA:chr01.1056");
		assertNotNull(p.getGoResults());
		assertEquals(2, p.getGoResults().size());
		assertTrue(p.getGoResults().contains("GO:0006355"));
		assertTrue(p.getGoResults().contains("GO:0043401"));
	}
}
