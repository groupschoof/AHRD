package ahrd.test;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

import ahrd.controller.Utils;

public class UtilsTest {

	@Test
	public void testFromFile() throws IOException {
		List<String> fromFile = Utils.fromFile("test/resources/ground_truth_go_annotations.goa");
		assertEquals(6, fromFile.size());

		fromFile = Utils.fromFile("./test/resources/empty_file.txt");
		assertEquals(0, fromFile.size());
	}

	@Test
	public void testRandomTrueOrFalse() {
		List<Boolean> rands = new ArrayList<Boolean>();
		for (int i = 0; i < 1000000; i++) {
			boolean r = Utils.randomTrueOrFalse();
			assertTrue("Random zero-or-one should be either TRUE or FALSE.",
					(r || !r));
			rands.add(r);
		}
		// Should contain at least a single 1 and a single zero:
		assertTrue("One in 1,000,000 random TrueOrFalse should be True",
				rands.contains(true));
		assertTrue("One in 1,000,000 random TrueOrFalse should be False",
				rands.contains(false));
	}

	@Test
	public void testRandomMultipleOfTen() {
		List<Long> rands = new ArrayList<Long>();
		for (int i = 0; i < 100000; i++) {
			Long r = Utils.randomMultipleOfTen();
			assertTrue(
					"Random Multiple-Of-Ten should fullfill >=10 and <=100, but is "
							+ r, (r >= 10 && r <= 100));
			rands.add(r);
		}
		Set<Long> distRands = new HashSet<Long>(rands);
		// Should have 10 entries:
		assertEquals(10, distRands.size());
		// Should contain one of each 10,20,30,..,100:
		for (long t = 10; t <= 100; t += 10) {
			assertTrue("Random Multiple-Of-Ten should contain " + t,
					distRands.contains(new Long(t)));
		}
	}

	@Test
	public void testRandomMultipleOfOneTenth() {
		List<Double> rands = new ArrayList<Double>();
		for (int i = 0; i < 100000; i++) {
			Double r = Utils.randomMultipleOfOneTenth();
			assertTrue(
					"Random Multiple-Of-One-Tenth should fullfill >=0.1 and <=1.0, but is "
							+ r, (r >= 0.1 && r <= 1.0));
			rands.add(r);
		}
		Set<Double> distRands = new HashSet<Double>(rands);
		// Should have 10 entries:
		assertEquals(10, distRands.size());
		// Should contain one of each 10,20,30,..,100:
		for (long t = 10; t <= 100; t += 10) {
			Double d = new Double(t * 0.01);
			assertTrue("Random Multiple-Of-One-Tenth should contain " + d,
					distRands.contains(d));
		}
	}

}
