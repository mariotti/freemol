!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H MODULE sections (Freemol F.Mariotti 03 Gen 2001)
!H-----------------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H 
!H This module manage the sections.
!H It is a part of the file manager modules system.
!H File manager moldule system is composed by:
!H basic I/O: baseio, Extended File manager: extio, Section Manager: sections.
!H 
!H This module define routines to handle sections
!H 
!H 
!H 
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H 
!
module sections
  !
  use vartypes
  use messages
  use extio
  !
  implicit none
  private
  !
!
!H
!H A proposal for functions in this module
!H 
!H  section_name	Get the name of the current section
!H  section_level	Get the level of the current section
!H  section_handler	Return an handler to the section
!H  
!H  
!H  
!H  section_next	make current the next section
!H  section_previous    make current the next section
!H  
!H  
!H  So propose again the question!!
!H  What do we do with a section??
!H  
!H  Seeking for data.
!H  Asking for format.
!H  
!H  
!H  
!H  
!H  The Section Handler
!H
!H This is the new type definition for the section handler  
!H
!H   type section
!H      ! This section
!H      character(FLCHARS) :: name
!H      integer(FINT) :: id
!H      integer(FINT) :: level
!H      integer(FINT) :: nrec
!H      integer(FINT) :: dim
!H      integer(FINT) :: numsub
!H      ! The file
!H      type(extfile) :: secfile
!H      ! The next section
!H      type(section), pointer :: nextsection
!H      character(FLCHARS) :: nextsectionname
!H      integer(FINT) :: nextsectionid
!H      ! The previous section
!H      type(section), pointer :: previoussection
!H      character(FLCHARS) :: previoussectionname
!H      integer(FINT) :: previoussectionid
!H      ! The section parameters
!H      character(FLCHARS), dimension(3,:) :: parameters
!H   end type section
!H  
!H  The data in this new type are:
!H  name	the name of this section
!H  id		the id of this section
!H  level	the level of this section
!H  nrec	the number of record for the section
!H              (WARNING: these can format dependent!)
!H  dim		This is the dimension of the section
!H              but is not yet implemented:
!H              It can be the number of lines in an ascii file
!H              or the file size in a binary file ....
!H  numsub	The number of subsections in this section
!H  secfile	A pointer to file handler used to determine
!H              the routines to be used to search inside the file
!H              and determine how to use the file.
!H  nextsection 	It is a pointer to the next section
!H  nextsectionname	the name of the next section
!H  nextsectionid 	the next section id
!H  previoussection 	It is a pointer to the previous section
!H  previoussectionname	the name of the previous section
!H  previoussectionid 	the previous section id
!H  subsection		It is a pointer to the subsections
!H			and actually the first subsection is
!H			the parameters subsection which is
!H			always builted unless we are parameters ;) !
!H  subsectionname	The name of the first subsection.
!H  subsectionid	the id of the subsection
!H  parameters          This vector contain the section parameters
!H                      it is a string vector of dimension (*,3).
!H                      The first column contain the parameter name
!H                      while the second column contain the parameter
!H                      value and the last one contain the type of the
!H                      parameter. A first idea about the parameter type
!H                      format can be:
!H                      X FD {K}
!H                      Where X is one of I,F,C meaning respectively
!H                      Integer, Floating point and Character and
!H                      FD is the dimension of the value filed and
!H                      optionally the last term is the KIND.
!H                      The actual form of the parameter type field
!H                      is not of particular interest because the
!H                      section module will implement methods to retrive
!H                      this data.
!H
!H  
!H The module can be actually used in two ways:
!H We can use always the calls with section handlers
!H or we can use the calls by names.
!H The main difference is: Calling by name we have the sections
!H uncorrelated but we can easily jump and skip
!H from section to section (this metod gives more freedom
!H but is supposed to be slower).
!H Using the structure with section handlers require
!H a bit more of work but it is supposed to be faster.
!H Actually the binary file can be handled only thrue the section
!H manager with this last method.
!H  
!H  
!H  
  type section
     ! This section
     character(FLCHARS) :: name
     integer(FINT) :: id
     integer(FINT) :: level
     integer(FINT) :: nrec
     integer(FINT) :: dim
     integer(FINT) :: numsub
     ! The file
     type(extfile) :: secfile
     ! The next section
     type(section), pointer :: nextsection
     character(FLCHARS) :: nextsectionname
     integer(FINT) :: nextsectionid
     ! The previous section
     type(section), pointer :: previoussection
     character(FLCHARS) :: previoussectionname
     integer(FINT) :: previoussectionid
     ! the subsections handler
     type(section), pointer :: subsection
     character(FLCHARS) :: subsectionname
     integer(FINT) :: subsectionid
  end type section
  !
  !
  integer(FINT), save, public :: sections_error = 0
  logical, save, private :: sections_isinit

  !
  ! Interfaces
  !-----------
  ! Mainly some sections procedure can be called
  ! by name (string) or by section(section type).
  !
contains
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H Function integer(FINT) function sections_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function sections_init()
    ! define the head section
    !
    sections_error = 0
    sections_isinit = .true.
    !
    sections_init = 0
    !
  end function sections_init
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H Function integer(FINT) function section_openfile(iun,filename,iotype,mainsec)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function section_openfile(iun,filename,iotype,mainsec)
    !
    integer(FINT), intent(out) :: iun
    character(FLCHARS), intent(in) :: filename
    type(section), intent(inout), optional :: mainsec
    character(FLCHARS), intent(in), optional :: iotype
!
!H  
!H   We make this supposition:
!H   If you call this procedure without the
!H   mainsec parameter you are going to open an ascii file.
!H   the other optional parameter iotype will define the use
!H   of the file.
!H  
!H   So the procedure will have to work like this:
!H   a call with only iun and filename require
!H   you are open an ascii file and it is in reading/append
!H   mode if the file exist otherwise it is opened for writing.
!H  
!H   if you call the procedure with only the extra parameter mainsec
!H   we require
!H   if the file exist it is opened in the current mode of the file
!H   (ascii or binary) but if the file does not exist it is
!H   opened as binary.
!H  
!H   the extra parameter iotype will override the default:
!H   the keywords and as you will see later the string to be
!H   given are as follow:
!H  
!H    R		The file is opened for reading
!H    W		the file is opened for writing
!H    U		Update or append at the end of the file
!H    N		New file: WARNING: it overwrite an existing file!
!H  
!H    A		Ascii file
!H    B		Binary file
!H    P		Portable file
!H    C		Ascii Compressed file
!H    T		naTive mode
!H    X		eXtra definiftion
!H  
!H   the strings are case insensitive and you can use the short
!H   version as well the long version, here some details:
!H  
!H   The given string of words is first checked for extended words
!H   than for two character words and finally for the short version
!H   which has to be one word string.
!H  
!H   if out of ambiguity mixed versions are allwed.
!H  
!H   Here is guide for matching:
!H  
!H   R = READ = RO(ReadOnly)
!H   W = WRITE = WO(WriteOnly which is actually append on existing files).
!H   U = UPDATE = UD(UpDate): in principle is supposed to update
!H                            the file and to create new sections as
!H                            needed but is not yet implemented.
!H   N = NEW = no two words version: The file is supposed to be or
!H                                   to became new!
!H  		                      which actually meas to overwrite
!H                                   the existing file.
!H   A = ASCII = no two words version: We force to read the file in ascii frmt!
!H   B = BINARY = no two words version: We force the file to be binary.
!H   P = PORTABLE = no two words version: It will be the default portable frmt!
!H   C = COMPRESSED = no two words version: It is the gzip compressed file.
!H                                          WARNING: it is possible we will
!H                                          use the bzip format instead!
!H   T = NATIVE = no two words version: it is the Freemol native file which
!H                                      will actually match a known onw
!H                                      but we leave alive this possiblity!
!H   X = EXTRA = no two words version: It has to be user/programmer defined
!H                                     WARNING: it is not yet implemented.
!H   The Algoritm will be:
!H  
!H   (NOTE: we are doing std io so we are not supposed to stress performance
!H          at least opening section files!)
!H  
!H   If it is only one word we check the word against the
!H   full version:
!H     if it is the full keyword we just use it with defaults.
!H     if it is not we perform interpretaion character by character
!H        in order to get the short version
!H   If it is more words it cannot be the compressed version
!H      so we check word by  word agaist (first the complete version
!H      and after against the two char version.
!H   For each operation we use the module extio!
!H   Keep in mind: sections files are not allowed to be SCRATCH
!H                 or anyway temporary files.
!H   For temporary storage of data we will implement special calls
!H   to procedures.
!H 
!
    character(FLCHARS) :: lfilename
    type(section) :: lmainsec
    character(FLCHARS) :: liotype
    !
    integer(FINT) :: irc
    !
    !Initialization
    section_openfile = 0
    iun=0
    lfilename = trim(filename)
    !
    ! We are asked for a file so
    ! we do it by extio module
    !-------------------------
    if(present(mainsec)) then
       lmainsec = mainsec ! I guess we can remove it ...
       irc = 0 !extio_open(mainsec%
    end if
    ! WARNING WE RETURN ERROR BECAUSE WE ARE NOT DONE
    section_openfile = -1
    
  end function section_openfile
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H Function integer(FINT) function section_next()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function section_next()
    !
    section_next = -1
    !
  end function section_next
  !
  !
  !
end module sections
