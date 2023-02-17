! Module file degree_trig.f90:
! From:
! http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi/0640/bks/SGI_Developer/books/MProF90_PG/sgi_html/apb.html#id35104

module degree_trig
   real(16), parameter :: &
      quadpi = 3.141592653589793238462643383279502884197Q0
   real(16), parameter :: dgr_to_rad = (quadpi/180Q0)
   intrinsic cos, sin, tan
contains
   function sind(dgr_argument)
   real(4) sind, dgr_argument
      sind = sin(dgr_to_rad * dgr_argument)
   end function

   function cosd(dgr_argument)
   real(4) cosd, dgr_argument
      cosd = cos(dgr_to_rad * dgr_argument)
   end function

   function tand(dgr_argument)
   real(4) tand, dgr_argument
      tand = tan(dgr_to_rad * dgr_argument)
   end function

   function dsind(dgr_argument)
   real(8) dsind, dgr_argument
      dsind = sin(dgr_to_rad * dgr_argument)
   end function

   function dcosd(dgr_argument)
   real(8) dcosd, dgr_argument
      dcosd = cos(dgr_to_rad * dgr_argument)
   end function

   function dtand(dgr_argument)
   real(8) dtand, dgr_argument
      dtand = tan(dgr_to_rad * dgr_argument)
   end function
end ! module
