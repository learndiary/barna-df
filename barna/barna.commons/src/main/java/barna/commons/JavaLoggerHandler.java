/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.commons;

import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;

/**
 * Delegates log messages to our {@link Log} implementations
 *
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */
public class JavaLoggerHandler extends Handler{
    @Override
    public void publish(LogRecord record) {
        Level level = record.getLevel();
        if(level.intValue() == Level.INFO.intValue()){
            Log.info(record.getLoggerName(), record.getMessage());
            return;
        }else if(level.intValue() == Level.WARNING.intValue()){
            Log.warn(record.getLoggerName(), record.getMessage());
        }else if(level.intValue() == Level.SEVERE.intValue()){
            Log.error(record.getMessage(), record.getThrown());
        }else {
            Log.debug(record.getLoggerName(), record.getMessage(), record.getThrown());
        }
    }

    @Override
    public void flush() {
    }

    @Override
    public void close() throws SecurityException {
    }
}
