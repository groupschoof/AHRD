package ahrd.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import ahrd.controller.Evaluator;
import ahrd.exception.MissingAccessionException;
import ahrd.exception.MissingProteinException;

public class EvaluatorTest {

	private Evaluator evaluator;

	@Before
	public void setup() throws IOException {
		this.evaluator = new Evaluator("./test/resources/ahrd_input.yml");
	}

	@Test
	public void testSetupBlast2GoAnnots() throws IOException,
			MissingAccessionException {
		evaluator.initializeProteins();
		evaluator.setupCompetitors();
		assertEquals(3, evaluator.getProteins().size());
		assertNotNull(
				"After setting up Competitors the Evaluator should have assigned a CompetitorAnnotation to the protein with accession 'gene:chr01.1056:mRNA:chr01.1056'.",
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
						.getEvaluationScoreCalculator().getCompetitorAnnotations());
		assertNotNull(
				"After setting up Blast2GoAnnots the Evaluator should have assigned a Blast2GoAnnot to the protein with accession 'gene:chr01.502:mRNA:chr01.502'.",
				evaluator.getProteins().get("gene:chr01.502:mRNA:chr01.502")
						.getEvaluationScoreCalculator().getCompetitorAnnotations());
		assertEquals("nrpb6a dna binding dna-directed rna polymerase",
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
						.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getDescription());
		assertEquals(6,
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
						.getEvaluationScoreCalculator().getCompetitorAnnotations().get("blast2go").getEvaluationTokens().size());
	}

	@Test
	public void testSetupReferenceDescriptions() throws IOException,
			MissingAccessionException {
		evaluator.initializeProteins();
		evaluator.setupReferenceDescriptions();
		assertEquals(3, evaluator.getProteins().size());
		assertNotNull(evaluator.getProteins()
				.get("gene:chr01.1056:mRNA:chr01.1056")
				.getEvaluationScoreCalculator().getReferenceDescription());
		assertNotNull(evaluator.getProteins()
				.get("gene:chr01.502:mRNA:chr01.502")
				.getEvaluationScoreCalculator().getReferenceDescription());
		assertEquals("Receptor-like protein kinase", evaluator.getProteins()
				.get("gene:chr01.1056:mRNA:chr01.1056")
				.getEvaluationScoreCalculator().getReferenceDescription()
				.getDescription());
		assertEquals(4,
				evaluator.getProteins().get("gene:chr01.1056:mRNA:chr01.1056")
						.getEvaluationScoreCalculator()
						.getReferenceDescription().getTokens().size());
	}

	@Test
	public void testTokenizesUnchangedBlastResults() throws IOException,
			MissingAccessionException, MissingProteinException, SAXException {
		evaluator.initializeProteins();
		evaluator.parseBlastResults();
		// TEST:
		assertEquals(
				"Parsing the BlastResults should also trigger tokenization of the Description-Lines of the unchanged BlastHits!",
				4, evaluator.getProteins().get("gene:chr01.502:mRNA:chr01.502")
						.getEvaluationScoreCalculator()
						.getUnchangedBlastResults().get("swissprot")
						.getTokens().size());
	}
	
	@Test
	public void testSetupCompetitors() throws IOException, MissingAccessionException {
		this.evaluator = new Evaluator("./test/resources/evaluator_example.yml");
		evaluator.initializeProteins();
		evaluator.setupReferenceDescriptions();
		evaluator.setupCompetitors();
	}
}
