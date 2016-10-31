# freemol tools manual

# Install

git clone https://github.com/mariotti/freemol.git

# First run: check install

Change to the freemol directory (_cd freemol/Freemol_) and type:

    ./config/autoconfigure

If everything is fine it prints out few lines like this:

    Got machine:  m_generic_osx
     Available compilers for m_generic_osx are:
    (  gfortran  )
    Testing Presence...
    Cheching for: gfortran
    WARNING: The autoconfigure script is not updated to local configuration
     Using gfortran as build compiler
    ./config/autoconfigure: line 167: [: ==: unary operator expected
    We should run as:
    ./config/configure m_generic_osx gfortran /Users/mariotti/GIT/freemol/Freemol

The first 3 lines depends on a predefined configuration. If you get an error
there you might need to create a configuration for your machine.

The script then checks for the presence of a pre-configured fortran compiler.

The last line tells you how you should run the configure program.

