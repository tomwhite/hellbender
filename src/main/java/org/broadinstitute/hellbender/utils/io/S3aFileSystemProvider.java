package org.broadinstitute.hellbender.utils.io;

public class S3aFileSystemProvider extends com.upplication.s3fs.S3FileSystemProvider {
    @Override
    public String getScheme() {
        return "s3a";
    }
}
