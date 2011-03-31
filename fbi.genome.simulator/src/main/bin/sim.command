scriptdir="`dirname ${0}`"
scriptdir=`cd "$scriptdir"; pwd`
java -Xmx1G -XX:-UseBiasedLocking -Xdock:name="Flux Simulator" -Xdock:icon="$scriptdir/../pics/favicon_mac.png" -DwrapperDir="$scriptdir" -jar "$scriptdir/../lib/FluxSimulator.jar" $@
