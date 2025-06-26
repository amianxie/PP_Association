# Protein-Protein Association Rate Constant Prediction 2015.pl
#open
#
use strict;
use LWP::Simple;

my (@list,@pro,@in,%noa1,%noa2,%nor1,%nor2,%dx1,%dx2,%dr1,%dr2,%ion,%xi,%ch1,%ch2,%exp,$nord)=();


my $home = "./";          #

my $conn_depend = 'd';
my $parameter = 0;
my $str_trj = 0;
my ($tfi, $sfi,$bfi, $pdbid,$c1r,$c2r,$c1l,$c2l,$r1,$r2,$Enat,$Esurf,$date, $ss1, $ss2) =();
print "\n";

if($ARGV[0])
{
    &process_parameter;
    #debug#print "process parameter "
}
else
{
    $nord = 10000;
    print "How many trajectories do you want to collect? (Default: 10000)  ";
    chomp ($nord=<>);
    print "PDB ID: (e.g. 1e4k) ";
    chomp ($list[0]=<>);
    print "Chian no. 1: (e.g. AB) ";
    chomp ($ch1{$list[0]}=<>); 
    print "Chian no. 2: (e.g. C) ";
    chomp ($ch2{$list[0]}=<>);
    print "How many residues are on the Chian $ch1{$list[0]}? ";
    chomp ($nor1{$list[0]}=<>); 
    print "How many residues are on the Chian $ch2{$list[0]}? ";
    chomp ($nor2{$list[0]}=<>); 
    print "How many atoms are on the Chian $ch1{$list[0]}? ";
    chomp ($noa1{$list[0]}=<>); 
    print "How many atoms are on the Chian $ch2{$list[0]}? ";
    chomp ($noa2{$list[0]}=<>);
    print "ionic strength (mM)? ";
    chomp ($ion{$list[0]}=<>);
    &mw2dxr;
}

my @addp; #addition properties
my $term = 'parc2015';


my $rehersal = 0;
my $repeat = 1;
my $cleaning = 0;
my $diff_xi = 1;
my $check = 1;
my $randomize_doing = 1;
my $large_scale = 0;
my $location = "everywhere"; #"ae/gordon"; #
my $broken_node = 0;
my $pause = 0;
my $emergent = 0;
my $jbnm = "";

my $H = 20;
my ($n1, $n2, $upp1, $upp2);

&showtime;

my $escape = 0;
my $overall_result = $home."/results/results_".$date.".txt";
open TOT, ">>$overall_result" or die "cannot open $overall_result:$!";
print TOT "$date\n";
foreach my $j(0..$nord-1)
{
    my $fno = $j;
    
    $addp[0] = "./parc2015";
    $addp[0] = "./parc2015trj" if $str_trj;
    print "Trajectory No. $j \n";
    #<>;
    
    @pro=@list;
    if ($randomize_doing )
    {
        my $noe = @pro;
        foreach my $i(0..@pro-1)
        {
            my $r=$i;
            if ($noe > 1) #while($r==$fno)
            {
                $r=int(rand($noe));
                redo if $r == $i;
                
                my $temp = $pro[$i];
                $pro[$i]=$pro[$r];
                $pro[$r]=$temp;
                
                $r=int(rand($noe));
                redo if $r == $i;
                
                $temp = $pro[$i];
                $pro[$i]=$pro[$r];
                $pro[$r]=$temp;
            }
        }
    }
    foreach my $i(0..@pro-1)
    {
        print TOT "$i $pro[$i]\n" unless $j; #debug#
        if ($j==0)
        {
            $tfi = $home.$pro[$i].'.pdb';
            $sfi = $home.$pro[$i]."_".$ch1{$pro[$j]}."_".$ch2{$pro[$j]}.'.pdb';
            if(!-e $sfi)
            {
                print "doesn't exist $sfi ";
                $parameter = $i;
                &desire2detail;
            }
            &first_round; #compute native & surf E
        }
        $addp[1] = $home.$pro[$i]."_".$ch1{$pro[$i]}."_".$ch2{$pro[$i]}.".pdb";        
        $addp[3] = $home."/trj/trj_".$pro[$i]."_".$fno.".pdb"; #$addp[3] = $home."trj/trj_".$pro[$i]."_".."_".$fno.".pdb";
        $addp[2] = $home."/rec/st_".$pro[$i]."_".$fno.".txt";  #$addp[2] = $home."stat/st_".$pro[$i]."_".."_".$fno.".txt";
        $addp[4] = $home."/ene/ene_".$pro[$i]."_".$fno.".txt"; #$addp[4] = $home."ene/ene_".$pro[$i]."_".$fno.".txt";
        $addp[5] = $noa1{$pro[$i]};
        $addp[6] = $noa2{$pro[$i]};
        $addp[7] = $nor1{$pro[$i]};
        $addp[8] = $nor2{$pro[$i]};
        $addp[9] = $dx1{$pro[$i]};
        $addp[10]= $dx2{$pro[$i]};
        $addp[11]= $dr1{$pro[$i]};
        $addp[12]= $dr2{$pro[$i]};
        $addp[15]= $xi{$pro[$i]} if $diff_xi; #xi
        unless(-s $addp[2] & -s $addp[4])
        {
            system "rm $addp[2] $addp[3] $addp[4]\n" if ($check and (-e $addp[2] or -e $addp[4]));
            print "@addp\n";
            system "@addp\n";
            die "The program didn't execute correctly. Please check the error message. If you cannot identify the error, please kindly forward all the error message to Paul Xie at amianxie\@gmail.com.\n" unless (-e $addp[2] & -e $addp[4]);
        }
    }
    
    #<>;
    &showtime;
}
print TOT "$nord trajectories have been recorded & collected.\n";
&parse_result;
close TOT;
#end of main




sub showtime
{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =localtime(time);
    $mon+=1;
    $year-=100;
    
    #print "TIME: $mday/$mon/$year,$isdst";
    $mon = "0".$mon if ($mon < 10);
    $mday= "0".$mday if ($mday<10);
    $date = $year.$mon.$mday;
    print "TIME: $mday/$mon/$year,$hour:$min:$sec\n";
}
sub process_parameter
{
    &default_value;
    print "default_value parameter$parameter @ARGV\n"; #debug#
    foreach my $i(0..@ARGV-1)
    {
        if($ARGV[$i]=~m/\-trj/)
        {
            $nord=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-pdbid/)
        {
            $list[$parameter]=$ARGV[$i+1];
            $ch1{$list[$parameter]} =$ARGV[$i+2];
            $ch2{$list[$parameter]} =$ARGV[$i+3];
        }
        if($ARGV[$i]=~m/\-w1/)
        {
            $addp[13]=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-w2/)
        {
            $addp[14]=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-T/)
        {
            $addp[16]=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-xi/)
        {
            $xi{$list[$i]}=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-atom1/)
        {
            $noa1{$list[$parameter]}=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-atom2/)
        {
            $noa2{$list[$parameter]}=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-res1/)
        {
            $nor1{$list[$parameter]}=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-res2/)
        {
            $nor2{$list[$parameter]}=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-dx1/)
        {
            $dx1{$list[$parameter]}=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-dx2/)
        {
            $dx2{$list[$parameter]}=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-dr1/)
        {
            $dr1{$list[$parameter]}=$ARGV[$i+1];
        }
        if($ARGV[$i]=~m/\-dr2/)
        {
            $dr2{$list[$parameter]}=$ARGV[$i+1];
        }
    }
    unless($dx1{$list[$parameter]}>0. and $dx2{$list[$parameter]}>0. and $dr1{$list[$parameter]}>0. and $dr2{$list[$parameter]}>0.)
    {
        $parameter = 0;
        &mw2dxr;
    }
}
sub default_value
{
    $nord=10000;
    $list[$parameter]='1brs';#'1e4k';
    $ch1{$list[$parameter]}='A';#'AB';
    $ch2{$list[$parameter]}='D';#'C';
    $addp[13] = 1.;  #we
    $addp[14] = .04;  #wh
    $addp[15]= 9.5; #7.6 #xi
    $addp[16]= 5.;  #temp
    $noa1{$list[$parameter]}=864; #3438
    $noa2{$list[$parameter]}=693; #1384
    $nor1{$list[$parameter]}=108; #432
    $nor2{$list[$parameter]}=87;  #172
    $dx1{$list[$parameter]}=0.;#11.7; #7.4
    $dx2{$list[$parameter]}=0.;#12.6; #10.0
    $dr1{$list[$parameter]}=0.;#4.6;  #1.1
    $dr2{$list[$parameter]}=0.;#5.7;  #2.9
    $xi{$list[$parameter]}=9.5;
}
sub first_round
{
    $addp[1] = $home.$pro[$parameter]."_".$ch1{$pro[$parameter]}."_".$ch2{$pro[$parameter]}.".pdb";     
    $addp[3] = "a3.pdb";#$home."natE/trj_".$pro[$parameter]."_natE.pdb"; #$addp[3] = $home."trj/trj_".$pro[$parameter]."_".$fno.".pdb";
    $addp[2] = "a2.txt";#$home."natE/st_".$pro[$parameter]."_natE.txt";  #$addp[2] = $home."stat/st_".$pro[$parameter]."_".$fno.".txt";
    $addp[4] = $home."natE/ene_".$pro[$parameter]."_natE.txt"; #$addp[4] = $home."ene/ene_".$pro[$parameter]."_".$fno.".txt";
    $addp[5] = $noa1{$pro[$parameter]};
    $addp[6] = $noa2{$pro[$parameter]};
    $addp[7] = $nor1{$pro[$parameter]};
    $addp[8] = $nor2{$pro[$parameter]};
    $addp[9] = $dx1{$pro[$parameter]};
    $addp[10]= $dx2{$pro[$parameter]};
    $addp[11]= $dr1{$pro[$parameter]};
    $addp[12]= $dr2{$pro[$parameter]};
    $addp[15]= $xi{$pro[$parameter]} if $diff_xi; #xi
    
    $addp[4] = $home."surfE/ene_".$pro[$parameter]."_surfE.txt"; #$addp[4] = $home."ene/ene_".$pro[$parameter]."_".$fno.".txt";
    if (!-e $addp[4])
    {
        system "./parc_surfE $addp[1] $addp[2] $addp[3] $addp[4] $addp[5] $addp[6] $addp[7] $addp[8] $addp[9] $addp[10] $addp[11] $addp[12] $addp[13] $addp[14] $addp[15] $addp[16]";
        print "./parc_surfE $addp[1] $addp[2] $addp[3] $addp[4] $addp[5] $addp[6] $addp[7] $addp[8] $addp[9] $addp[10] $addp[11] $addp[12] $addp[13] $addp[14] $addp[15] $addp[16]\n";
        #<>;
    }
    open IN, "<$addp[4]" or die "cannot open $addp[4]:$!";
    my $sen=<IN>;
    close IN;
    if($sen =~m/0 (\-?\d+\.\d+) /)
    {#0 -43.1165 -43.1165 0 nan 90 8.19131e-15 -12 -17.2 -42.4285 0 1
        $Esurf = $1;
        print "electrostatic potential of the binding surface: $Esurf\n";
        print TOT "electrostatic potential of the binding surface: $Esurf\n";
    }
    else
    {
        print "data error! Please check $addp[4]!\n";
    }
    
    $addp[4] = $home."natE/ene_".$pro[$parameter]."_natE.txt"; #$addp[4] = $home."ene/ene_".$pro[$parameter]."_".$fno.".txt";
    if (!-e $addp[4])
    {
        system "./parc_natE $addp[1] $addp[2] $addp[3] $addp[4] $addp[5] $addp[6] $addp[7] $addp[8] $addp[9] $addp[10] $addp[11] $addp[12] $addp[13] $addp[14] $addp[15] $addp[16]";
        print "./parc_natE $addp[1] $addp[2] $addp[3] $addp[4] $addp[5] $addp[6] $addp[7] $addp[8] $addp[9] $addp[10] $addp[11] $addp[12] $addp[13] $addp[14] $addp[15] $addp[16]\n";
        #<>;
    }
    open IN, "<$addp[4]" or die "cannot open $addp[4]:$!";
    $sen=<IN>;
    close IN;
    if($sen =~m/0 (\-?\d+\.\d+) /)
    {#0 -43.1165 -43.1165 0 nan 90 8.19131e-15 -12 -17.2 -42.4285 0 1
        $Enat = $1;
        print "electrostatic potential of the whole protein: $Enat\n";
        printf "the ratio r(elec): %7.4f\n",$Esurf/$Enat;
        print TOT "electrostatic potential of the whole protein: $Enat\n";
        printf TOT "the ratio r(elec): %7.4f\n",$Esurf/$Enat;
    }
    else
    {
        print "data error! Please check $addp[4]!\n";
    }
    
    my $output3 = $home."SS/".$pro[$parameter]."_".$ch1{$pro[$parameter]}."_".$ch2{$pro[$parameter]}."_SSc.txt";
    if (!-e $output3)
    {
        open SS, ">$output3" or die "cannot qopen $output3 :$!";# if(!-e $output2);
        #foreach my $i(0..@list-1)
        {
            $parameter = $pro[$parameter];
            print "node3: $parameter"; #debug
            &ssratio;
            #<>;
        }
        close SS;
        #die; #debug#
    }
    open IN, "<$output3" or die "cannot open $output3:$!";
    $sen=<IN>;
    close IN;
    if($sen =~m/\w{4} \d+ \d+ \d+ \d+ (\d+\.\d+) \d+ \d+ \d+ \d+ (\d+\.\d+)/)
    {#1brs 108 26 16 10 38.4615 87 24 11 13 54.1667 0 "%d %d %d %d %7.4f %d\n"
        ($ss1, $ss2) = ($1,$2);
        print "%interface residues on the flexible loops of chain $ch1{$pro[$parameter]}: $ss1 %\n";
        print "%interface residues on the flexible loops of chain $ch2{$pro[$parameter]}: $ss2 %\n";
        print TOT "%interface residues on the flexible loops of chain $ch1{$pro[$parameter]}: $ss1 %\n";
        print TOT "%interface residues on the flexible loops of chain $ch2{$pro[$parameter]}: $ss2 %\n";
    }
    else
    {
        print "data error! Please check $output3! If you cannot identify the error, please kindly forward this error message and the file $output3 to Paul Xie at amianxie\@gmail.com.\n";
    }
    
    
}#end first_round
sub ssratio
{
    my %charge=('ASP',1,'GLU',1,'ARG',1,'LYS',1,'HIS',1);
    my $pdb = $home.$parameter."_".$ch1{$parameter}."_".$ch2{$parameter}.".pdb";
    `ls -l $pdb`;
    #next;
    
    my (@tbs, @tbe, @ths, @the, @tts, @tte, @cont)=();
    my ($ntbs,$ntbe,$nths,$nthe,$ntts,$ntte,)=(0,0,0,0,0,0);
    my (@mbs, @mbe, @mhs, @mhe, @mts, @mte)=();
    my ($nmbs,$nmbe,$nmhs,$nmhe,$nmts,$nmte,)=(0,0,0,0,0,0);
    my $content ='';
    my $url = '';
    
    my $loc_chain = 1;
    my ($thr, $dist_thr)=(0,6);
    my $desire_chains1 = length $ch1{$parameter};
    my $desire_chains2 = length $ch2{$parameter};
    print "SS: chain details: $parameter, $ch1{$parameter} $desire_chains1 $ch2{$parameter} $desire_chains2\n";
    
    foreach my $j(1..$desire_chains1)
    {
        $loc_chain = 1;
        my @eachchain = split//,$ch1{$parameter};
        my $c101;
        
        while ($loc_chain)
        {
            #4debug#print "pro_part1: $desire_chains1 chain(s) $ch1{$parameter} \n";
            $url = 'http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode='.$parameter.'&template=protein.html&o=SUMMARY&l='.$loc_chain.'&c=1&chain='.$ch1{$parameter};
                    #'http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode='.$parameter.'&template=protein.html&o=STRANDS&l=2&s=1&c=6&chain=A
            #4debug#print $url,"\n";
            
            #4debug#<>;
            $content = get($url);
            print "$content\n";
            die "Couldn't get $url" unless defined $content;
            @cont = split(/\<\/a\>\<BR\>/,$content);
            #4debug#print $content;
            $c101='';
            if ($cont[0]=~m/\<IMG align\=bottom width\=18 border\=0 src=\"\/thornton\-srv\/databases\/pdbsum\/templates\/gif\/chain($eachchain[$j-1])\.jpg\"\>/)#[\s^(]{0,}\((\d+) residues\)#($cont[$i]=~m/chain\_on\.src=\"\/thornton\-srv\/databases\/PDBsum.*$parameter\/\/chain(\w)/)
            {#Chain <IMG align=bottom width=18 border=0 src="/thornton-srv/databases/pdbsum/templates/gif/chainA.jpg"> &
             #Chain <IMG align=bottom width=18 border=0 src="/thornton-srv/databases/pdbsum/templates/gif/chainL.jpg"> (33 residues)#chain_on.src="/thornton-srv/databases/PDBsum/jw/1jwh//chainA.jpg";
                #4debug#print "note2 $loc_chain, $j, $eachchain[$j-1], $1,\n";
                $c101 = $1;
                last;
                #<>;
            }
            elsif ($cont[0]=~m/\<IMG align\=bottom width\=18 border\=0 src=\"\/thornton\-srv\/databases\/pdbsum\/templates\/gif\/chain(\w?)\.jpg\"\>/)
            {
                $c101 = $1;
                if ($ch1{$parameter} eq $c101)
                {
                    print "chain$loc_chain $c101 eq chain1 $ch1{$parameter}!\n";
                    
                }
                elsif ($ch1{$parameter}=~m/$c101/)
                {
                    print "chain$loc_chain $c101 =~ chain1 $ch1{$parameter}!\n";
                    #<>;
                }
                elsif ($ch2{$parameter}=~m/$c101/)
                {
                    
                    print "but chain$loc_chain $c101 =~ chain2 $ch2{$parameter}!! \n";
                    #unshift @cont, $c101;
                    #next;
                }
                else
                {
                    print "weird! please check: chain$loc_chain $c101 ne chain1 $ch1{$parameter} or chain2 $ch2{$parameter}!!!\n";
                    <>;
                }
                
            }
            else
            {
                die "cannot find the information of Secondary structure of $parameter on PDBsum: $url\n";
            }
            $loc_chain++;
            $loc_chain=1 if $loc_chain > $desire_chains1;
        }
        foreach my $i (0..20) #(0..@cont-1)
        {
            #4debug#print "j:$j,i$i, ";
            if($cont[$i]=~m/\<\/body\>/)
            {
                last;
            }
            
            if ($cont[$i]=~m/(\d+) strand/)
            {
                #4debug#print $cont[$i],$i,"\n";
                #4debug#print "chain1 $i, $1 strand(s)\n";
                #4debug#<>;
                my @cell = split(/\<\/td\>/,$cont[$i+1]);
                
                foreach my $j (1..@cell-1)
                {
                    #4debug#print "$cell[$j],$j,note!";
                    #4debug#<>;
                    if($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+[A-Z]?)\<\/font\>\&nbsp\;/)#($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+)[ABCDEFGHIJKLMN]?\<\/font\>\&nbsp\;/)
                    {
                        if($ntbs>$ntbe)
                        {
                            $tbe[$ntbe]=$1.$c101;
                            $ntbe++;
                            #4debug#print "end $ntbe $tbe[$ntbe-1]\n";
                            #4debug#<>;
                        }
                        else #ntbs nof T beta start
                        {
                            $tbs[$ntbs]=$1.$c101;
                            $ntbs++;
                            #4debug#print "start $ntbs $tbs[$ntbs-1]\n";
                            #4debug#<>;
                        }
                    }
                }
                #4debug#<>;
            }
            elsif ($cont[$i]=~m/(\d+) helice/ or $cont[$i]=~m/(1) helix/)
            {
                #4debug#print $cont[$i],$i,"\n";
                #4debug#print "chain1 $i, $1 helice(s)\n";
                #4debug#<>;
                my @cell = split(/\<\/td\>/,$cont[$i+1]);
                
                foreach my $j (1..@cell-1)
                {
                    #4debug#print "$cell[$j],$j,note!";
                    #4debug#<>;
                    if($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+[A-Z]?)\<\/font\>\&nbsp\;/)#($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+)[ABCDEFGHIJKLMN]?\<\/font\>\&nbsp\;/)
                    {
                        if($nths>$nthe)
                        {
                            $the[$nthe]=$1.$c101;
                            $nthe++;
                            #4debug#print "end $nthe $the[$nthe-1]\n";
                            #4debug#<>;
                        }
                        else #nths nof T helix start
                        {
                            $ths[$nths]=$1.$c101;
                            $nths++;
                            #4debug#print "start $nths $ths[$nths-1]\n";
                            #4debug#<>;
                        }
                    }
                }
                #4debug#<>;
            }
            elsif ($cont[$i]=~m/(\d+) beta turn/)
            {
                #4debug#print $cont[$i],$i,"\n";
                #4debug#print "chain1 $i, $1 beta turn(s)\n";
                #4debug#<>;
                my @cell = split(/\<\/td\>/,$cont[$i+1]);
                
                foreach my $j (1..@cell-1)
                {
                    #4debug#print "$cell[$j],$j,note!";
                    #4debug#<>;
                    if($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+[A-Z]?)\-([A-Z][a-z][a-z]\d+[A-Z]?)\<\/font\>\&nbsp\;/)#($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+)[ABCDEFGHIJKLMN]?\-([A-Z][a-z][a-z]\d+)[ABCDEFGHIJKLMN]?\<\/font\>\&nbsp\;/)
                    {
                        $tts[$ntts]=$1.$c101; #resnm
                        $ntts++; #ntts nof T bturn start
                        #4debug#print "start $ntts $tts[$ntts-1] ";
                        $tte[$ntte]=$2.$c101;
                        $ntte++;
                        #4debug#print "end $ntte $tte[$ntte-1]\n";
                        #4debug#<>;
                    }
                }
                #4debug#<>;
            }
        }
        #4debug#<>;
    }
    #4debug#print "\n";
    $loc_chain=0 if $loc_chain != $desire_chains1; #$loc_chain=1 if $loc_chain != $desire_chains1;
    foreach my $j(1..$desire_chains2)
    {
        $loc_chain++;#$loc_chain = 1;
        my @eachchain = split//,$ch2{$parameter};
        my $c102;
        
        while ($loc_chain)
        {
            #4debug#print "pro_part2: $desire_chains2 chain(s) $ch2{$parameter} \n";
            $url = 'http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode='.$parameter.'&template=protein.html&o=SUMMARY&l='.$loc_chain.'&c=1&chain='.$ch2{$parameter};
                     #'http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode='.$parameter.'&template=protein.html&o=STRANDS&l=2&s=1&c=6&chain=A
            #4debug#print $url,"\n";
            #use LWP::Simple;
            #<>;
            $content = get($url);
            die "Couldn't get $url" unless defined $content;
            @cont = split(/\<\/a\>\<BR\>/,$content);
            #4debug#print $content;
            $c102='';
            if ($cont[0]=~m/\<IMG align\=bottom width\=18 border\=0 src=\"\/thornton\-srv\/databases\/pdbsum\/templates\/gif\/chain($eachchain[$j-1])\.jpg\"\>/)#[\s^(]{0,}\((\d+) residues\)#($cont[$i]=~m/chain\_on\.src=\"\/thornton\-srv\/databases\/PDBsum.*$parameter\/\/chain(\w)/)
            {#Chain <IMG align=bottom width=18 border=0 src="/thornton-srv/databases/pdbsum/templates/gif/chainA.jpg"> &
             #Chain <IMG align=bottom width=18 border=0 src="/thornton-srv/databases/pdbsum/templates/gif/chainL.jpg"> (33 residues)#chain_on.src="/thornton-srv/databases/PDBsum/jw/1jwh//chainA.jpg";
                #4debug#print "note1 $loc_chain, $j, $eachchain[$j-1], $1,\n";
                $c102 = $1;
                last;
                #4debug#<>;
            }
            elsif ($cont[0]=~m/\<IMG align\=bottom width\=18 border\=0 src=\"\/thornton\-srv\/databases\/pdbsum\/templates\/gif\/chain(\w?)\.jpg\"\>/)
            {
                $c102 = $1;
                $loc_chain=0 if ($c102 eq '');
            }
            else
            {
                die "cannot find the information of Secondary structure of $parameter on PDBsum: $url\n";
            }
            $loc_chain++;
            my $tmp=$desire_chains1+$desire_chains1;
            $loc_chain=1 if ($loc_chain > $tmp);#$desire_chains12);
        }
        
        foreach my $i (0..20) #(0..@cont-1)
        {
            #4debug#print "j$j,i(2)$i, ";
            if($cont[$i]=~m/\<\/body\>/)
            {
                last;
            }
            
            if ($cont[$i]=~m/(\d+) strand/)
            {
                #4debug#print $cont[$i],$i,"\n";
                #4debug#print "chain2 $i, $1 strand(s)\n";
                #4debug#<>;
                my @cell = split(/\<\/td\>/,$cont[$i+1]);
                
                foreach my $j (1..@cell-1)
                {
                    #4debug#print "$cell[$j],$j,note!";
                    #4debug#<>;
                    if($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+[A-Z]?)\<\/font\>\&nbsp\;/)#($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+)[ABCDEFGHIJKLMN]?\<\/font\>\&nbsp\;/)
                    {
                        if($nmbs>$nmbe)
                        {
                            $mbe[$nmbe]=$1.$c102;
                            $nmbe++;
                            #4debug#print "end $nmbe $mbe[$nmbe-1]\n";
                            #4debug#<>;
                        }
                        else
                        {
                            $mbs[$nmbs]=$1.$c102;
                            $nmbs++;
                            #4debug#print "start $nmbs $mbs[$nmbs-1]\n";
                            #4debug#<>;
                        }
                    }
                }
                #4debug#<>;
            }
            elsif ($cont[$i]=~m/(\d+) helice/ or $cont[$i]=~m/(1) helix/)
            {
                #4debug#print $cont[$i],$i,"\n";
                #4debug#print "chain2 $i, $1 helices\n";
                #4debug#<>;
                my @cell = split(/\<\/td\>/,$cont[$i+1]);
                
                foreach my $j (1..@cell-1)
                {
                    #4debug#print "$cell[$j],$j,note!";
                    #4debug#<>;
                    if($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+[A-Z]?)\<\/font\>\&nbsp\;/)#($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+)[ABCDEFGHIJKLMN]?\<\/font\>\&nbsp\;/)
                    {
                        if($nmhs>$nmhe)
                        {
                            $mhe[$nmhe]=$1.$c102;
                            $nmhe++;
                            #4debug#print "end $nmhe $mhe[$nmhe-1]\n";
                            #4debug#<>;
                        }
                        else
                        {
                            $mhs[$nmhs]=$1.$c102;
                            $nmhs++;
                            #4debug#print "start $nmhs $mhs[$nmhs-1]\n";
                            #4debug#<>;
                        }
                    }
                }
                #4debug#<>;
            }
            elsif ($cont[$i]=~m/\d+ beta turn/)
            {
                #4debug#print $cont[$i],$i,"\n";
                #4debug#print "chain2 $i, $1 beta turn(s)\n";
                #4debug#<>;
                my @cell = split(/\<\/td\>/,$cont[$i+1]);
                
                foreach my $j (1..@cell-1)
                {
                    #4debug#print "$cell[$j],$j,note!";
                    #4debug#<>;
                    if($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+[A-Z]?)\-([A-Z][a-z][a-z]\d+[A-Z]?)\<\/font\>\&nbsp\;/)#($cell[$j]=~m/\<font face=arial\>([A-Z][a-z][a-z]\d+)[ABCDEFGHIJKLMN]?\-([A-Z][a-z][a-z]\d+)[ABCDEFGHIJKLMN]?\<\/font\>\&nbsp\;/)
                    {
                        $mts[$nmts]=$1.$c102;
                        $nmts++;
                        #4debug#print "start $nmts $mts[$nmts-1] ";
                        $mte[$nmte]=$2.$c102;
                        $nmte++;
                        #4debug#print "end $nmte $mte[$nmte-1]\n";
                        #4debug#<>;
                    }
                }
                #4debug#<>;
            }
        }
    }
    #4debug#print "nodeSS2: $parameter\n";
    #4debug#<>;
    
    my (@tss,@mss,%tindx,%mindx,%taa,%maa)=();
    
    my $infile = $home.$parameter."_".$ch1{$parameter}."_".$ch2{$parameter}.".pdb";  #"./".$parameter.".pdb";
    open IN, "<$infile" or die "cannot open $infile:$!";
    my @in=<IN>;
    close IN;
    print "open $infile: \n";
    my (@trid,@tx,@ty,@tz,@tsurf,@trnm,@tbeta,@mrid,@mx,@my,@mz,@msurf,@mrnm,@mbeta) = ();
    my ($mhc,$tcr)=();
    my ($mna,$tna,$mstart,$mstop,$tstart,$tstop,$rindx,$mem_c,$mem_r)=(0,0,0,0,0,0,0,'',''); 
    foreach my $i(0..@in-1)
    {
        if($in[$i]=~m/ATOM\s{1,6}\d{1,6} [\w\s]\w{1,3}\s{0,2}[\s1A]?(\w{3}) (\w)\s{0,3}((-?\d{1,4})[A-Z]?)\s{3,4}(\s{0,3}-?\d{1,4}\.\d{3})(\s{0,3}-?\d{1,4}\.\d{3})(\s{0,3}-?\d{1,4}\.\d{3})  \d\.\d\d(\s{0,2}\d{1,3}\.\d{2})/)
        {#most typical 20150512
         #ATOM      1  N   VAL A   3      16.783  48.812  26.447  1.00 30.15           N
         #ATOM   5654  N   GLU C   6      -9.387  11.979 -62.127  1.00110.01           N
         #ATOM    832  NH1AARG A 113      82.049  42.594  18.825  0.50  2.00           N
         #ATOM   4348  N   HIS B  -1      61.255 -26.547   4.711  1.00 68.77           N
            my ($aa,$chain,$resid,$resno,$x,$y,$z,$b)=($1,$2,$1.$3,$4,$5,$6,$7,$8); #resid
            $resid.=$chain;
            if ($ch1{$parameter}=~m/$chain/)#($chain eq $ch1{$parameter})
            {
                if($i==0)#($chain ne $mem_c)
                {
                    $mem_c = $chain;
                    $rindx = $resno;
                    #<>;
                }
                if($resid ne $mem_r)
                {
                    $tindx{$resid}=$rindx; #rindex
                    print "$aa,$chain,$rindx,$tindx{$resid},$resid($mem_r),$resno,$x,$y,$z\n";
                    $mem_r = $resid;
                    $rindx++;
                    #<>;
                }
                
                $tstart=$resno if $i == 0 ;
                ($tx[$tna],$ty[$tna],$tz[$tna],$trid[$tna],$tstop,$trnm[$tna],$tbeta[$tna])=($x,$y,$z,$resid,$resno,$aa,$b);
                $tna++; #tna, nof T atom
                #<>;
            }
            elsif($ch2{$parameter}=~m/$chain/)
            {
                if($ch1{$parameter}=~m/$mem_c/) #($chain ne $mem_c)
                {
                    $mem_c = $chain;
                    $rindx = $resno;
                }
                if($resid ne $mem_r)
                {
                    $mindx{$resid}=$rindx;
                    #print "$aa,$chain,$rindx,$mindx{$resid},$resid($mem_r),$resno,$x,$y,$z\n";
                    $mem_r = $resid;
                    $rindx++;
                    #<>;
                }
                
                $mstart=$resno if $mstart==0;
                ($mx[$mna],$my[$mna],$mz[$mna],$mrid[$mna],$mstop,$mrnm[$mna],$mbeta[$mna])=($x,$y,$z,$resid,$resno,$aa,$b);
                $mna++;
                #<>;
            }
            else
            {
                die "error2: $in[$i] $chain $ch1{$parameter} $ch2{$parameter}\n";
            }
        }
        elsif($in[$i]=~m/ATOM\s{1,6}\d{1,6} [\w\s]\w{1,3}\s{1,2}[\s23456789BCDEFGH]?(\w{3}) (\w)\s{1,3}((\d{1,4})[A-Z]?)\s{3,4}(\s{0,3}-?\d{1,4}\.\d{3})(\s{0,3}-?\d{1,4}\.\d{3})(\s{0,3}-?\d{1,4}\.\d{3})  \d\.\d\d(\s{0,2}\d{1,3}\.\d{2})/)
        {
        }
        elsif($in[$i]=~m/ATOM/)
        {
            die "error470: $in[$i]";
        }
        
    }
    {#2nd str
        
        foreach my $i(0..$ntbs-1)
        {
            my ($a,$b)=(-1,-1);
            $tbs[$i]=uc($tbs[$i]);
            $tbe[$i]=uc($tbe[$i]);
            $a=$tindx{$tbs[$i]};
            $b=$tindx{$tbe[$i]};
            if ($a ne '' and $b ne '')
            {
                print "BS $i $a $b $tbs[$i] $tbe[$i]\n";
                foreach my $j($a..$b)
                {
                    $tss[$j]="S";
                }
            }
        }
        foreach my $i(0..$nths-1)
        {
            my ($a,$b)=(-1,-1);
            $ths[$i]=uc($ths[$i]);
            $the[$i]=uc($the[$i]);
            $a=$tindx{$ths[$i]};
            $b=$tindx{$the[$i]};
            if ($a ne '' and $b ne '')
            {
                print "HX $i $a $b $ths[$i] $the[$i]\n";
                foreach my $j($a..$b)
                {
                    $tss[$j]="H";
                }
            }
        }
        foreach my $i(0..$ntts-1)
        {
            my ($a,$b)=(-1,-1);
            $tts[$i]=uc($tts[$i]);
            $tte[$i]=uc($tte[$i]);
            $a=$tindx{$tts[$i]};
            $b=$tindx{$tte[$i]};    
            if ($a ne '' and $b ne '')
            {
                print "BT $i $a $b $tts[$i] $tte[$i]\n";
                foreach my $j($a..$b)
                {
                    $tss[$j]="T";
                }
            }
        }
        
        foreach my $i(0..$nmbs-1)
        {
            my ($a,$b)=(-1,-1);
            $mbs[$i]=uc($mbs[$i]);
            $mbe[$i]=uc($mbe[$i]);
            $a=$mindx{$mbs[$i]};
            $b=$mindx{$mbe[$i]};
            if ($a ne '' and $b ne '')
            {
                print "BS $i $a $b $mbs[$i] $mbe[$i]\n";
                foreach my $j($a..$b)
                {
                    $mss[$j]="S";
                }
            }
        }
        foreach my $i(0..$nmhs-1)
        {
            my ($a,$b)=(-1,-1);
            $mhs[$i]=uc($mhs[$i]);
            $mhe[$i]=uc($mhe[$i]);
            $a=$mindx{$mhs[$i]};
            $b=$mindx{$mhe[$i]};
            if ($a ne '' and $b ne '')
            {
                print "HX $i $a $b $mhs[$i] $mhe[$i]\n";
                foreach my $j($a..$b)
                {
                    $mss[$j]="H";
                }
            }
        }
        foreach my $i(0..$nmts-1)
        {
            my ($a,$b)=(-1,-1);
            $mts[$i]=uc($mts[$i]);
            $mte[$i]=uc($mte[$i]);
            $a=$mindx{$mts[$i]};
            $b=$mindx{$mte[$i]};
            if ($a ne '' and $b ne '')
            {
                print "BT $i $a $b $mts[$i] $mte[$i]\n";
                foreach my $j($a..$b)
                {
                    $mss[$j]="T";
                }
            }
        }
    }        
    #<>;
    my $clash = 0;
    foreach my $i (0..$tna-1)
    {
        foreach my $j (0..$mna-1)
        {
            #next if ($tsurf[$tindx{$trid[$i]}] and $msurf[$mindx{$mrid[$j]}]);
            my $dist=(($tx[$i]-$mx[$j])**2+($ty[$i]-$my[$j])**2+($tz[$i]-$mz[$j])**2)**.5;
            if ($dist < 2)
            {
                $clash ++;
                print "clash!! $clash\n";
                <>;
            }
            if ($dist < $dist_thr) #contact surf
            {
                ($tsurf[$tindx{$trid[$i]}],$msurf[$mindx{$mrid[$j]}]) = (1,1);
                #print "surf $i, $j; $trid[$i], $mrid[$j]; $tindx{$trid[$i]}, $mindx{$mrid[$j]}; $tsurf[$tindx{$trid[$i]}], $msurf[$mindx{$mrid[$j]}]; $tss[$tindx{$trid[$i]}], $mss[$mindx{$mrid[$j]}]\n";
                #<>;
            }
        }
    }
    print "   12345678901234567890123456789012345678901234567890\n"; 
    foreach my $i($tstart/50..($tstart+$nor1{$parameter})/50)#(0..$n1{$parameter}/50)
    {
        printf "%3d",$i*50+1;
        foreach my $j(1..50)
        {
            $tss[$i*50+$j] = " " unless $tss[$i*50+$j];
            printf "%s",$tss[$i*50+$j];
        }
        print "\n";
        printf "%3d",$i*50+1;
        foreach my $j(1..50)
        {
            $tsurf[$i*50+$j] = " " unless $tsurf[$i*50+$j];
            printf "%s",$tsurf[$i*50+$j];
        }
        print "\n";
    }
    print "   12345678901234567890123456789012345678901234567890\n"; 
    foreach my $i($mstart/50..($mstart+$nor2{$parameter})/50)#(0..$n2{$parameter}/50)
    {
        printf "%3d",$i*50+1;
        foreach my $j(1..50)
        {
            $mss[$i*50+$j] = " " unless $mss[$i*50+$j];
            printf "%s",$mss[$i*50+$j];
        }
        print "\n";
        printf "%3d",$i*50+1;
        foreach my $j(1..50)
        {
            $msurf[$i*50+$j] = " " unless $msurf[$i*50+$j];
            printf "%s",$msurf[$i*50+$j];
        }
        print "\n";
    }
    
    my ($len,$tot,$tsurf,$msurf,$tcore,$mcore,$sh,$ss,$st,$sl,$ch,$cs,$ct,$cl,$beta,$sob,$cob)=($nor1{$parameter},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,);
    print "$tstart..($tstart+$nor1{$parameter}-1)\n";
    my ($prev,$tmpt,$tmpl,$tmpss) = (0,0,0,0);
    foreach my $i($tstart..($tstart+$nor1{$parameter}-1))#(1..$len)
    {
        #print "$ch1{$parameter}, $i, $tss[$i], $tsurf[$i], $tsurf,\n";
        #<>;
        $tot++;
        if($tsurf[$i] > 0) #surf or not
        {
            $tsurf++;
            #next unless $charge{$trnm[$i]};
            #next if (($i-$tstart)<1 or ($tstart+$nor1{$parameter}-$i) <2); #avoid N- C- Term
            if($tss[$i] eq 'H' or $tss[$i] eq 'S')
            {
                $tmpss++;
                #if($tmpl>0)
                #{
                    $prev=$tmpl+$tmpt;
                #}
                #else
                #{
                #    $prev=$tmpl;
                #}
                #the diff scenario
                #$sl = $prev if $prev > $sl;  
                $sl+= $tmpl+$tmpt; #if $prev > $thr;#if $tmpl > $thr;
                print " Loop (length):$sl\n" if $prev > 0;
                ($prev,$tmpl,$tmpt) = (0,0,0,);
            }
            else
            {
                if($tss[$i] eq ' ')
                {
                    $tmpl++;
                    print "$tss[$i]$i:";
                    #<>;
                }
                elsif ($tss[$i] eq 'T')
                {
                    $tmpt++;
                    print "$tss[$i]$i:";
                    #<>;
                }
                $ss+= $tmpss if $tmpss > $thr;
                print " SS (length):$ss\n" if $tmpss > $thr;
                ($tmpss) = (0,0,0,);
            }
            
            $beta+=$tbeta[$i];
            $sob+=$tbeta[$i];
        }
        else
        {
            #if($tmpl>0)
            #{
                $prev=$tmpl+$tmpt;
            #}
            #else
            #{
            #    $prev=$tmpl;
            #}
            #the diff scenario
            #$sl = $prev if $prev > $sl;
            $sl+= $tmpl+$tmpt;# if $prev > $thr;#if $tmpl > $thr;
            $ss+= $tmpss if $tmpss > $thr;
            print " loop (length):$sl\n" if $prev > 0;
            print " SS (length):$ss\n" if $tmpss > $thr;
            ($prev,$tmpl,$tmpt,$tmpss) = (0,0,0,0);
        }
    }
    print "$parameter $infile\n";
    #printf    "$parameter $tstart..$tstop %d %d %d %6.2f %6.2f %6.2f %6.2f %d %6.2f %6.2f %6.2f %6.2f \n",$len,$tot,$surf,$sh,$ss,$st,$sl,$core,$ch,$cs,$ct,$cl,;
    #<>;
    printf    "$parameter %d %d %d %d %d %7.4f \n",$len,$tot,$tsurf,$sl,$ss,$ss/$tsurf*100;
    printf SS "$parameter %d %d %d %d %7.4f ",$len,$tsurf,$sl,$ss,$sl/$tsurf*100,;
    #printf  "$parameter %d %d %d %d %d %d %d %d %d %d %d %d \n"  ,$len,$tot,$tsurf,$sh,$ss,$st,$sl,$tcore,$ch,$cs,$ct,$cl,;
    #printf  "$parameter %d %d %d %d %d %d %d %d %d %d %d %d "  ,$len,$tot,$tsurf,$sh,$ss,$st,$sl,$tcore,$ch,$cs,$ct,$cl,;
   
    ($len,$tot,$msurf,$mcore,$sh,$ss,$st,$sl,$ch,$cs,$ct,$cl,$prev, $tmpt, $tmpl,$tmpss)=($nor2{$parameter},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,);
    print "$mstart..($mstart+$nor2{$parameter}-1)\n";
    foreach my $i($mstart..($mstart+$nor2{$parameter}-1))#(1..$len)
    {
        #print "$ch2{$parameter}, $i, $mss[$i], $msurf[$i], $msurf,\n";
        $tot++;
        if($msurf[$i] > 0)
        {
            $msurf++;
            #next unless $charge{$mrnm[$i]};
            next if (($i-$mstart)<1 or ($mstart+$nor2{$parameter}-$i) <2);
            if($mss[$i] eq 'H' or $mss[$i] eq 'S')
            {
                $tmpss++;
                #if($tmpl>0)
                #{
                    $prev=$tmpl+$tmpt;
                #}
                #else
                #{
                #    $prev=$tmpl;
                #}
                #the diff scenario
                #$sl = $prev if $prev > $sl;
                $sl+= $tmpl+$tmpt;# if $prev > $thr;#if $tmpl > $thr;
                print " Loop (length):$sl\n" if $prev > 0;
                ($prev,$tmpl,$tmpt) = (0,0,0,);
            }
            else
            {
                if($mss[$i] eq ' ')
                {
                    $tmpl++;
                    print "$mss[$i]$i:";
                    #<>;
                }
                elsif ($mss[$i] eq 'T')
                {
                    $tmpt++;
                    print "$mss[$i]$i:";
                    #<>;
                }
                $ss+= $tmpss if $tmpss > $thr;
                print " SS (length):$ss\n" if $tmpss > $thr;
                ($tmpss) = (0,0,0,);
            }
            $beta+=$tbeta[$i];
            $sob+=$tbeta[$i];
        }
        else
        {
            #if($tmpl>0)
            #{
                $prev=$tmpl+$tmpt;
            #}
            #else
            #{
            #    $prev=$tmpl;
            #}
            #the diff scenario
            #$sl = $prev if $prev > $sl;
            $sl+= $tmpl+$tmpt;# if $prev > $thr;#if $tmpl > $thr;
            $ss+= $tmpss if $tmpss > $thr;
            print " loop (length):$sl\n" if $prev > 0;
            print " SS (length):$ss\n" if $tmpss > $thr;
            ($prev,$tmpl,$tmpt,$tmpss) = (0,0,0,0);
        }
    }
    #printf    "$parameter $mstart..$mstop %d %d %d %6.2f %6.2f %6.2f %6.2f %d %6.2f %6.2f %6.2f %6.2f \n",$len,$tot,$surf,$sh,$ss,$st,$sl,$core,$ch,$cs,$ct,$cl,;
    #<>;
    printf    "%d %d %d %d %d %7.4f clash $clash\n",$len,$tot,$msurf,$sl,$ss,$ss/$msurf*100,;
    printf SS "%d %d %d %d %7.4f %d\n",$len,$msurf,$sl,$ss,$sl/$msurf*100,;
    #printf SS "%d %d %d %d %d %d %d %d %d %d %d %d %6.2f %6.2f %6.2f\n",      $len,$tot,$msurf,$sh,$ss,$st,$sl,$mcore,$ch,$cs,$ct,$cl,$sob/($tsurf+$msurf),$cob/($tcore+$mcore),$beta/($tsurf+$msurf+$tcore+$mcore);
    #<>;
    #4debug#<>;
}#end SSratio
sub desire2detail
{   
    my $a = uc($pro[$parameter]);
    my  $url = 'http://www.rcsb.org/pdb/files/'.$a.'.pdb';
    my  $lfi = $home.$a.'.pdb';
        $tfi = $home.$pro[$parameter].'.pdb';
        $sfi = $home.$pro[$parameter]."_".$ch1{$pro[$parameter]}."_".$ch2{$pro[$parameter]}.'.pdb';
    if (!-e $tfi)
    {
        print "from $url 2 $lfi\n";
        getstore($url, $lfi);
        `mv $lfi $tfi`;
        ($c1l,$c2l,$r1,$r2)=(0,0,0,0);
        &transfer;
        $nor1{$pro[$parameter]}=$c1l;
        $nor2{$pro[$parameter]}=$c2l;
        $noa1{$pro[$parameter]}=$r1;
        $noa1{$pro[$parameter]}=$r2;
        
        #1KAC            73000           7.6     19.9    13.6    9.8     11.1    2.7     3.9     185     124     1401    959 #$pnm,$kon,$xi,$mw1,$mw2,$dx1,$dx2,$dr1,$dr2,$nor1,$nor2,$noa1,$noa2
        #print $i," $a,$pro[$parameter],$ch1{$pro[$parameter]},$ch2{$pro[$parameter]},$exp{$pro[$parameter]},$ion{$pro[$parameter]}\n";
        #<>;
    }
    
}
sub mw2dxr
{
    my ($mw1,$mw2);#dx=26.647/x^(1/3) dr=30.232/x*180/100
    
    $mw1=$nor1{$list[$parameter]}*.11;
    $mw2=$nor2{$list[$parameter]}*.11;
    $dx1{$list[$parameter]}=26.647/$mw1**(1/3);
    $dr1{$list[$parameter]}=30.232/$mw1*1.8/3.14;
    $dx2{$list[$parameter]}=26.647/$mw2**(1/3);
    $dr2{$list[$parameter]}=30.232/$mw2*1.8/3.14;
    $xi{$list[$parameter]}=(9241.6/$ion{$list[$parameter]})**.5 if ($ion{$list[$parameter]}>0.);# and $xi{$list[$parameter]} <= 0. );
    print "mw2dxr param:$parameter, pro:$list[$parameter], $mw1 $mw2 $dx1{$list[$parameter]} $dx2{$list[$parameter]} $dr1{$list[$parameter]} $dr2{$list[$parameter]} \n";
    #debug#<>;        
}
sub transfer
{
    print "transfer $tfi 2 $sfi\n";
    
    open IN, "<$tfi" or die "cannot open $tfi:$!";
    my @tfi = <IN>;
    close IN;
    
    my $parameter = 0;
    my @tmp = split(//,$ch1{$pdbid});
    $c1r = $tmp[0];
       @tmp = split(//,$ch2{$pdbid});
    $c2r = $tmp[0];
           
    #unless (-e $ccp and -e $cpr and -e $fp3 and -e $slg)
    {
        open PRO, ">$sfi" or die "cannot open $sfi: $!"; # unless -e $sfi;
        #open BK,  ">$bfi" or die "cannot open $bfi: $!"; # unless -e $sfi;
        
        my ($ter,$prev) = (0,-100);
        my %chain = ();
        my $no_chain = 0;
        
        
        foreach my $line (0..@tfi-1)
        {
            #print $line,$tfi[$line];
            if ($tfi[$line]=~m/^ATOM/)
            {
                #         1         2         3         4         5         6         7         8
                #12345678901234567890123456789012345678901234567890123456789012345678901234567890
                #HETATM 3186 ZN    ZN     1      30.355  20.927   2.323  1.00 50.77          ZN
                #ATOM   2292  N   ASN B  -1      29.225  25.915  -5.505  1.00104.99           N
                #ATOM     91  AD1 ASN A  14      34.700   7.470  18.600  1.00  0.00           .  
                #ATOM     31 1HE2 GLN L   3       9.560  -8.933  35.031  0.00  0.00           H
                #ATOM    137  N   SER A1000       9.760  35.069   7.428  1.00 29.03           N
                #ATOM      1  N   MET     1      24.461  59.918   4.546  1.00 28.24      5DFR N
                #ATOM    122  CG1AILE A  16      29.494  14.605  34.873  0.50 13.63           C
                #ATOM   1039 HH11BARG A  79      33.478  48.539  18.769  0.45 14.55           H
                #ATOM      2  C5'   C Q   0      12.941  12.992   5.572  1.00 75.84           C
                #ATOM   2312  P     A B   1      13.774  39.026  15.251  0.50220.71           P
                #ATOM   2315  O5'   A B   1      14.872  37.934  15.708  0.50218.41
                #ATOM   2336 H5''   A B   1      14.877  38.343  17.663  0.50210.20           H           O
                #ATOM   1506  CG AHIS A 501      13.993  16.566  16.039 -0.50  7.18           C
                #1 -  6        Record name   "ATOM  "
                #7 - 11        Integer       serial       Atom  serial number.
                #13 - 16        Atom          name         Atom name.
                #17             Character     altLoc       Alternate location indicator.
                #18 - 20        Residue name  resName      Residue name.
                #22             Character     chainID      Chain identifier.
                #23 - 26        Integer       resSeq       Residue sequence number.
                #27             AChar         iCode        Code for insertion of residues.
                #31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
                #39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
                #47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
                #55 - 60        Real(6.2)     occupancy    Occupancy.
                #61 - 66        Real(6.2)     tempFactor   Temperature  factor.
                #77 - 78        LString(2)    element      Element symbol, right-justified.
                #79 - 80        LString(2)    charge       Charge  on the atom.
                if   ($tfi[$line] =~ m/ATOM\s{2}\s{0,4}(\d{1,5})\s[\s1-4H](([A-Z\d]{1})[\w\s\'\"\*]{2})[BCDEFGH23456789]([\s\w]{3})\s([\s\w])\s{0,3}(\-?\d{1,4})[\w\s]\s{3}\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{2}(\d\.\d{2})\s{0,2}(\-?\d{1,3}\.\d{2})/)
                {#ATOM    188 1HG1BILE B  17      12.865  -2.614  29.863  0.74  7.56           H
                 #ATOM   1039 HH11BARG A  79      33.478  48.539  18.769  0.45 14.55           H
                    #$abmode{"$pdbid:$7:$8:$9"} ++;
                    #temp#print "AB model $pdbid $2, $4, $pdbid:$7:$8:$9 $tfi[$line]";#temp#
                    
                }
                elsif($tfi[$line] =~ m/ATOM\s{2}\s{0,4}(\d{1,5})\s([HD\d]([A-Z0-9])[\d\w\s\'\"\*]{2})[\sA1]([\s\w]{3})\s([\s\w])\s{0,3}(\-?\d{1,4})[\w\s]\s{3}\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{2}(\d\.\d{2})\s{0,2}(\-?\d{1,3}\.\d{2})/)
                {#ATOM    191 HG11 VAL A  12      -6.388   6.846  19.544  1.00 37.47
                    #temp#print "HD0 $pdbid $2, $4, $pdbid:$7:$8:$9 $tfi[$line]";#temp#
                    #<>;
                }
                elsif($tfi[$line] =~ m/ATOM\s{2}\s{0,4}(\d{1,5})\s\s(([HD])[\d\w\s\'\"\*]{2})[\sA1]([\s\w]{3})\s([\s\w])\s{0,3}(\-?\d{1,4})[\w\s]\s{3}\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{2}(\d\.\d{2})\s{0,2}(\-?\d{1,3}\.\d{2})/)
                {#ATOM    191 HG11 VAL A  12      -6.388   6.846  19.544  1.00 37.47
                    #temp#print "HD1 $pdbid $2, $4, $pdbid:$7:$8:$9 $tfi[$line]";#temp#
                    #<>;
                }
                elsif($tfi[$line] =~ m/ATOM\s{2}\s{0,4}\d{1,5}\s\s([123][HD][\d\w\s\'\"\*])[\sA1]([\s\w]{3})\s([\s\w])\s{0,3}(\-?\d{1,4})[\w\s]\s{3}\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{2}(\d\.\d{2})\s{0,2}(\-?\d{1,3}\.\d{2})\s{11}[HDZ]/)
                {#ATOM    697  1DZ LYS A  42      13.379  29.208  -7.309  1.00 24.10           H
                    #temp#print "HD0 $pdbid $2, $4, $pdbid:$7:$8:$9 $tfi[$line]";#temp#
                    #<>;
                }
                elsif($tfi[$line] =~ m/(ATOM\s{2}\s{0,4}(\d{1,5})\s\s(([A-Z\d])[\w\s\'\"\*]{2})[\sA1]([\w]{3})\s)([\w])(\s{0,3}(\-?\d{1,4})([\w\s\'])\s{3}\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})[\s\-]{0,2}(\d{1,2}\.\d{2})\s{0,2}(\-?\d{1,3}\.\d{2})(.*)$)/)
                {#                                    $2 serail   $3at nom$4element                $5 residue   $6chain id     $8 residue no                  $7 x                      $8 y                      $9 z                       $10 occup                $11 b-fac          
                    #ATOM    915  CB  LYS A 123      53.862  13.613   8.552-99.00 22.48           C
                    #ATOM      1  O5*   A C 901     -20.110  -9.664  -9.517  1.00 68.27           O
                    #ATOM   2437  N   LEU I   2'      3.580  36.459  12.382  1.00  0.19           N
                    my ($sen1,$sen7,$sr,$an,$ele,$res,$ci,$rno,$rno1,$x,$y,$z,$occ,$bf)=($1,$7,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,$13,$14);
                    my $rno0 =$rno.$rno1;
                    #print "$1$6$7\n$tfi[$line]";
                    #<>;
                    
                    if ($ch1{$pdbid} =~m/$ci/)
                    {
                        $parameter = 'reverse' if $ter > 1;
                        $ter = 1;
                        if($parameter eq 'reverse')
                        {
                            $r2++;
                        }
                        else
                        {
                            $r1++;
                        }
                        if($rno0 ne $prev)
                        {
                            if($parameter eq 'reverse')
                            {
                                $c2l++;
                            }
                            else
                            {
                                $c1l++;
                            }
                            $prev = $rno0;
                        }
                        #print PRO $sen1,$c1r,$sen7,"\n";
                        print PRO $tfi[$line];
                    }
                    elsif ($ch2{$pdbid} =~m/$ci/)
                    {
                        #next unless $ter > 0;
                        $parameter = 'reverse' if $ter < 1;
                        $ter = 2;
                        if($parameter eq 'reverse')
                        {
                            $r1++;
                        }
                        else
                        {
                            $r2++;
                        }
                        if($rno0 ne $prev)
                        {
                            if($parameter eq 'reverse')
                            {
                                $c1l++;
                            }
                            else
                            {
                                $c2l++;
                            }
                            $prev = $rno0;
                        }
                        #print PRO $sen1,$c2r,$sen7,"\n";
                        print PRO $tfi[$line];
                    }
                    #print PRO "$tfi[$line]";# if ($ter<3 and $ter >0);
                    
                    if ($ele eq 'H' or $ele eq 'D')
                    {
                        die "die-err01:$pdbid $tfi[$line]: $!";# if $ele eq 'H' or $ele eq 'D'; #temp#
                        return 0;
                        last;
                    }
                    
                }
                elsif($tfi[$line] =~ m/ATOM\s{2}\s{0,4}(\d{1,5})\s\s(([A-Z\d])[\w\s\'\"\*]{2})[\sA1]([\s\w]{3})\s([\s\w])\s{0,3}(\-?\d{1,4})([\w\s\'])\s{3}\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})\s{0,3}(\-?\d{1,4}\.\d{3})[\s\-]{0,2}(\d{1,2}\.\d{2})\s{0,2}(\-?\d{1,3}\.\d{2})(.*)/)
                {#
                    print "?$tfi[$line] ";
                    #<>;
                }
                else
                {
                    print "die-err02:$pdbid no match: $tfi[$line]\n";
                    return 0;
                    last;
                }
            }
            elsif ($tfi[$line]=~m/^HETATM/)
            {
            }
            elsif ($tfi[$line]=~m/^TER/)
            {
                
                if ($ter == 1)
                {
                    #print PRO "$tfi[$line]";
                    $prev = -100;
                }
                elsif ($ter == 2)
                {
                    #print PRO "$tfi[$line]";
                    $prev = -100;
                    #last;
                }
            }
            elsif ($tfi[$line]=~m/^ENDMDL/)
            {
                last;
            }
        }
        
        print "end transfer: $ter, $c1l, $c2l,\n";
        #<>;
    }
}
sub parse_result
{
    #my $ofi = "./r_".$term.".txt";
    #open OUT, ">$ofi" or die "cannot open $ofi:$!";
    
    foreach my $i(0..@list-1)
    {
        my ($model,$conn_depend,$thrr,$step,$max,$data_complete,$nofc,$detail,$thr0,$thr1,$thr2,$thr3,$avg2,$avg3,$avg4,$avg5,$avg6,$rho)=($parameter,'d',10,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,);
        my (@count,@t_array) = ();
        $model = $list[$i];
        print "($model) processing $nord trajectories ...";
        foreach my $j(0..$nord-1)
        {
            print ".";
            $data_complete = 0;
            my $data1 = $home."/trj/trj_".$model."_".$j.".pdb";
            my $data2 = $home."/rec/st_".$model."_".$j.".txt";
            my $data3 = $home."/ene/ene_".$model."_".$j.".txt";
            
            my $data7 = "test".$model;
            my ($suc,$nol2,$nol3,$nol4,$nol5)=(0,0,0,0,0);
            
            #print "node1: $i $data2 $data3 $data4 $data5 $data6 \n $nol2 $nol3 $nol4 $nol5 $avg2 $avg3 $avg4 $avg5 $avg6 \n";
            #<>;
            
            if($conn_depend eq 'd' or $conn_depend eq 'nr')
            {
                ($thr0,$thr1,$thr2,$thr3) = (-1,-2,-3,-7);
            }
            elsif($conn_depend eq 'eH' or $conn_depend eq 'n')
            {
                ($thr0,$thr1,$thr2,$thr3) = (-3,-4,-5,-6);
            }
            
            if (!-e $data3)
            {
                print "cannot find $data3!";
                next;
            }
            elsif (!-s $data3)
            {
                print "$data3 is empty!";
                next;
            }
            
            open IN, "<$data3" or die "cannot open $data3:$!";
            my @in = <IN>;
            close IN;
            my $des = $max/$step;
            $in[1] = $in[0] unless ($in[1]=~m/.{1,}/);
            if ($in[$des]=~m/$max\s\-?\d+\.?\d{0,}e?[+-]?\d{0,2}\s\-?\d+\.?\d{0,}e?[+-]?\d{0,2}\s\-?\d+\.?\d{0,}e?[+-]?\d{0,2}\s\-?\d+\.?\d{0,}e?[+-]?\d{0,2}\s\-?\d+\.?\d{0,}e?[+-]?\d{0,2}\s(\-?\d+\.?\d{0,}e?[+-]?\d{0,2})\s(-?\d+)\s(\-?\d+\.?\d{0,}e?[+-]?\d{0,2})\s(\-?\d+\.?\d{0,}e?[+-]?\d{0,2})\s(\-?\d+\.?\d{0,}e?[+-]?\d{0,2})\s/) #($in[$desire]=~m/\d\s9999\s\d+\.\d+\s(-?\d+)\s/)
            {# 9999 -0.00542038 -0.00542038 54.9339 94.5931 93.6646 18.8745 0 0 -0.00542038 0
                #1000 -6.08565e-06 -6.08565e-06 101.731 118.217 82.8148 46.737 0 0 -6.08565e-06 0
             #  1 693 12.2327 -1 -1.8 -36.1665 0
                $data_complete++;
                my ($rmsd,$t_ass,$en) = ($1,$2,$4);
                {
                    if ($t_ass < $thr0 and $rmsd < $thrr) #-1
                    {
                        $count[1]++;
                        if ($t_ass < $thr1 and $rmsd < $thrr) #-2
                        {
                            $count[2]++; #this is the answer
                            if ($t_ass < $thr2 and $rmsd < $thrr) #-3
                            {
                                $count[3]++;
                                if ($t_ass < $thr3 and $rmsd < $thrr) #-7
                                {
                                    $count[4]++;
                                }
                            }
                            if($detail)
                            {
                                print "asso: $j,  $rmsd, $t_ass, $en, $in[$des]";
                                <>;
                            }
                        }
                    }
                    if($detail and $rmsd < 10 and $en < -10)
                    {
                        print "$model, $j,  $rmsd, $t_ass, $en, $in[$des]";
                        <>;
                    }
                }
                #printf "$i $max $desire %d, %d, %d, %d, %d, $in[$desire] ", $t_ass, $count[1], $count[2],$count[3], $count[4];
                #<>;
                #$avg = ($avg*($count-1)+$t_ass)/$count;
                #$t_array[$j]=$t_ass;
                #print "$i.";
                #print "node2: $i $count[0] $count[1] $t_ass $avg2,$avg3,$avg4,$avg5,$avg6 \n";
            }
            else
            {
                print "error1 $data3 $j $in[$des] \n";
                system "rm $data1 $data2 $data3" if $cleaning;
                next;
                die;
                <>;
            }
            open IN, "<$data2" or die "cannot open $data2:$!";
            @in = <IN>;
            close IN;
            if ($in[0]=~m/(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s\-?\d+\.?\d{0,}e?[+-]?\d{0,2}\s(-?\d+)\s/) #($in[$desire]=~m/\d\s9999\s\d+\.\d+\s(-?\d+)\s/)
            {#1 1 0 0 0 0 0 10001 9364 9825 9359 2.48074 -12 -8.4 -19.4238 0
            #ofi<<nof_trj<<" "<<mem[10]<<" "<<mem[0]<<" "<<mem[1]<<" "<<mem[2]<<" "<<mem[3]<<" "<<mem[4]<<" "<<mem[11]<<" "<<mem[12]<<" "<<mem[13]<<" "<<mem[14]<<" "<<trj[0]<<" "<<trj[1]<< " " <<trj[2]<< " " <<trj[3]<< " " <<trj[4]<< endl;
            #10:ever re-associate > 3, 0~4: last stat > 0~4, 11~14: total steps;move fw steps;total steps bf re-ass; move fw stps bf re-ass
                ($suc,$nol2,$nol3,$nol4,$nol5)= ($2,$8,$9,$10,$11);
                $data_complete++;
            }
            else
            {
                print "error2 $data2 $j $in[0] \n";
                system "rm $data1 $data2 $data3" if $cleaning;
                next;
                die "error2  $data2 $j $in[0] ";
                <>;
            }
            $nofc++ if $data_complete > 1;
            #$avg2=($avg2*$j+$nol2)/($j+1);
            #$avg3=($avg3*$j+$nol3)/($j+1);
            #$avg4=($avg4*$j+$nol4)/($j+1);
            #$avg5=($avg5*$j+$nol5)/($j+1);
            #$avg6=($avg6*)/($i+1);
            #$count[11]++ if $nol4 < $max;
            #$t_array[$j]=$nol4;
        }
        
        my @new = sort {$b <=> $a} @t_array;
        #foreach my $i(0..$desire-1)
        #{
        #    print "$new[$i] ";
        #}
        #<>;
        #print "\n";
        $rho = $count[2]/$nofc;
        printf "\nPDB id: $model; Max. duration: $max ns; no. trj: $nofc/$nord; No. of successfully associated: $count[2]\n";#%3d %3d %3d %3d avg %5.2f %5.2f %5.2f %5.2f %5.2f medium %5.1f max %5d \n",$count[11], $count[1], $count[2],$count[3], $count[4], $avg2, $avg3, $avg4, $avg5, $avg6, ($new[499]+$new[500])/2, $new[0];
        printf "Predicted kon of $model: %e /Ms\n", $rho/.00167/(1-$rho)/$max/.000000001;
        printf TOT "PDB id: $model; Max. duration: $max ns; no. trj: $nofc/$nord; No. of successfully associated: $count[2]\n";#%3d %3d %3d %3d avg %5.2f %5.2f %5.2f %5.2f %5.2f medium %5.1f max %5d \n",$count[11], $count[1], $count[2],$count[3], $count[4], $avg2, $avg3, $avg4, $avg5, $avg6, ($new[499]+$new[500])/2, $new[0];
        printf TOT "Predicted kon of $model: %e /Ms\n", $rho/.00167/(1-$rho)/$max/.000000001;
    }
    #close OUT;
}
