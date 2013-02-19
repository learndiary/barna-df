package barna.commons.parameters;

import org.junit.Test;

import java.util.EnumSet;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * Created with IntelliJ IDEA.
 * User: Emilio Palumbo
 * Date: 6/19/12
 * Time: 6:05 PM
 */
public class EnumSetParameterTest {

    private enum A{A,B,C}

    @Test
    public void testEnumParamCreation() throws Exception {
        EnumSetParameter<A> p = new EnumSetParameter<A>("My param", "",EnumSet.noneOf(A.class),A.class,null);
        p.set(EnumSet.noneOf(A.class));
        p.parse("B");
        assertEquals(1, p.get().size());
        assertTrue(p.get().contains(A.B));
    }
}
