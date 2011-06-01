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

package fbi.commons;

/**
 * @deprecated will bre removed soon, do not use !
 */
public interface ReadyOrNot {

    public boolean isReady();

    public void set(Object settings);

    public boolean loadStats();

    public boolean isFinished();

    public void setLoadStats(boolean val);

    public boolean getLoadStats();

    public void killResult();
}
