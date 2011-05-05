package fbi.commons.parameters;

import java.io.File;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface FileNameParser {
    File parse(String string);
}
