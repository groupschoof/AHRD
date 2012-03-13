package ahrd.test;

import static org.junit.Assert.*;
import org.junit.Test;
import ahrd.controller.AHRD;

public class AhrdTestRun {

	public AhrdTestRun() {
		super();
	}

	@Test
	public void testRun() throws Exception {
		String[] args = {"./test/resources/ahrd_input_test_run.yml"};
		AHRD.main(args);
		// No error?, then we are happy..
		assertTrue(true);
	}
}
