#!/usr/bin/perl

use 5.10.1 ;
#opens file and transfers it's contents to $string.
open $outfile, './out' or die "I couldn't get at out.txt: $!";
open (SURVIVAL, '>>survivalrates.txt') ;

$string = <$outfile> ;

@lineparts = split( / / , $string ) ;
$value = @lineparts[3] ;
print SURVIVAL $value ;

while($string = <$outfile>)
{

	if($string =~ /Selected popn/i )
	{
		@lineparts = split( / / , $string ) ;
		$value = @lineparts[3] ;
		print SURVIVAL $value ;
	}
}

	close $outfile ;
        close (SURVIVAL) ;




