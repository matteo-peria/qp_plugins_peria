module io_color

  implicit none

  !> Single precision real numbers
  integer, parameter :: sp = selected_real_kind(6)

  !> Double precision real numbers
  integer, parameter :: dp = selected_real_kind(15)

  !> Char length for integers
  integer, parameter :: i1 = selected_int_kind(2)

  !> Short length for integers
  integer, parameter :: i2 = selected_int_kind(4)

  !> Length of default integers
  integer, parameter :: i4 = selected_int_kind(9)

  !> Long length for integers
  integer, parameter :: i8 = selected_int_kind(18)


  !> Container for terminal escape code
  type :: color_code
    !> Style descriptor
    integer(i1) :: style = -1_i1
    !> Background color descriptor
    integer(i1) :: bg = -1_i1
    !> Foreground color descriptor
    integer(i1) :: fg = -1_i1
  end type color_code

  !> Colorizer class for handling colorful output in the terminal
  type, public :: color_output

    type(color_code) :: &
      reset = color_code(), &
      bold = color_code(), &
      dim = color_code(), &
      italic = color_code(), &
      underline = color_code(), &
      blink = color_code(), &
      reverse = color_code(), &
      hidden = color_code()

    type(color_code) :: &
      black = color_code(), &
      red = color_code(), &
      green = color_code(), &
      yellow = color_code(), &
      blue = color_code(), &
      magenta = color_code(), &
      cyan = color_code(), &
      white = color_code()

    type(color_code) :: &
      bg_black = color_code(), &
      bg_red = color_code(), &
      bg_green = color_code(), &
      bg_yellow = color_code(), &
      bg_blue = color_code(), &
      bg_magenta = color_code(), &
      bg_cyan = color_code(), &
      bg_white = color_code()
  end type color_output

  interface color_output
    module procedure :: new_color_output
  end interface color_output

  ! Modified from private to public
  type(color_output), public, protected :: color


contains



!> Transform a color code into an actual ANSI escape sequence
pure function escape_color(code) result(str)
  !> Color code to be used
  type(color_code), intent(in) :: code
  !> ANSI escape sequence representing the color code
  character(len=:), allocatable :: str
  character, parameter :: chars(0:9) = &
    ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]

  if (anycolor(code)) then
    str = achar(27) // "[0"  ! Always reset the style
    if (code%style > 0 .and. code%style < 10) str = str // ";" // chars(code%style)
    if (code%fg >= 0 .and. code%fg < 10) str = str // ";3" // chars(code%fg)
    if (code%bg >= 0 .and. code%bg < 10) str = str // ";4" // chars(code%bg)
    str = str // "m"
  else
    str = ""
  end if
end function escape_color


!> Check whether the code describes any color or is just a stub
pure function anycolor(code)
  !> Escape sequence
  type(color_code), intent(in) :: code
  !> Any color / style is active
  logical :: anycolor

  anycolor = code%fg >= 0 .or. code%bg >= 0 .or. code%style >= 0
end function anycolor

  !> Initialize color output
  subroutine init_color_output(use_color)
    !> Enable color output
    logical, intent(in) :: use_color

    color = new_color_output(use_color)
  end subroutine init_color_output

  !> Create a new colorizer object
  function new_color_output(use_color) result(new)
    !> Enable color output
    logical, intent(in) :: use_color
    !> New instance of the colorizer
    type(color_output) :: new

    type(color_code), parameter :: &
      reset = color_code(style=0_i1), &
      bold = color_code(style=1_i1), &
      dim = color_code(style=2_i1), &
      italic = color_code(style=3_i1), &
      underline = color_code(style=4_i1), &
      blink = color_code(style=5_i1), &
      reverse = color_code(style=7_i1), &
      hidden = color_code(style=8_i1)

    type(color_code), parameter :: &
      black = color_code(fg=0_i1), &
      red = color_code(fg=1_i1), &
      green = color_code(fg=2_i1), &
      yellow = color_code(fg=3_i1), &
      blue = color_code(fg=4_i1), &
      magenta = color_code(fg=5_i1), &
      cyan = color_code(fg=6_i1), &
      white = color_code(fg=7_i1)

    type(color_code), parameter :: &
      bg_black = color_code(bg=0_i1), &
      bg_red = color_code(bg=1_i1), &
      bg_green = color_code(bg=2_i1), &
      bg_yellow = color_code(bg=3_i1), &
      bg_blue = color_code(bg=4_i1), &
      bg_magenta = color_code(bg=5_i1), &
      bg_cyan = color_code(bg=6_i1), &
      bg_white = color_code(bg=7_i1)

    if (use_color) then
      new%reset = reset
      new%bold = bold
      new%dim = dim
      new%italic = italic
      new%underline = underline
      new%blink = blink
      new%reverse = reverse
      new%hidden = hidden
      new%black = black
      new%red = red
      new%green = green
      new%yellow = yellow
      new%blue = blue
      new%magenta = magenta
      new%cyan = cyan
      new%white = white
      new%bg_black = bg_black
      new%bg_red = bg_red
      new%bg_green = bg_green
      new%bg_yellow = bg_yellow
      new%bg_blue = bg_blue
      new%bg_magenta = bg_magenta
      new%bg_cyan = bg_cyan
      new%bg_white = bg_white
    end if
  end function new_color_output


end module io_color
