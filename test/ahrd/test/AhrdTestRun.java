package ahrd.test;

import static org.junit.Assert.*;
import org.junit.Test;
import ahrd.controller.AHRD;
import ahrd.controller.Evaluator;

public class AhrdTestRun {

	public AhrdTestRun() {
		super();
	}

	@Test
	public void testRun() throws Exception {
		String[] args = { "./test/resources/ahrd_input_test_run.yml" };
		AHRD.main(args);
		// No error?, then we are happy..
		assertTrue(true);
	}

	@Test
	public void testRunWithDomainArchitectureSimilarityScoring()
			throws Exception {
		String[] args = { "./test/resources/ahrd_evaluator_dom_arch_sim_input.yml" };
		Evaluator.main(args);
		// No error?, then we are happy..
		assertTrue(true);
	}
}
