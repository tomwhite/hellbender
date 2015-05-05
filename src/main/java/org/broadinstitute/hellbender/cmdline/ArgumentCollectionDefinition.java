package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.Serializable;

/**
 * marker interface for classes that are intended to be used with @ArgumentCollection
 * Those are parsed by the CommandLineParser class.
 */
public interface ArgumentCollectionDefinition extends Serializable{

    /**
     * Implementing classes can an override this in order to provide custom argument validation that is more complicated
     * than what can be enforced by the argument parser.
     * @throws UserException.BadArgumentValue
     */
    default void validate() throws UserException.BadArgumentValue{
        //defaults to valid
    }

}
