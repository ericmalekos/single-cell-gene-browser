-- Launcher for the single-cell gene browser. install_mac.sh fills in the __PLACEHOLDERS__ and compiles this
-- into "Gene Browser.app".
-- Double-click: start the server if it is not running, then open the app in the default browser.
-- Double-click while it runs: open it again, or stop the server.

property scbrowserPath : "__SCBROWSER__"
property serverPort : "__PORT__"
property logPath : "__LOG__"
property pidPath : "__PIDFILE__"

on appURL()
	return "http://localhost:" & serverPort & "/app"
end appURL

on serverUp()
	try
		do shell script "/usr/bin/curl -s -o /dev/null -m 2 " & quoted form of (my appURL())
		return true
	on error
		return false
	end try
end serverUp

on startServer()
	-- the redirections let do shell script return at once while the server keeps running
	do shell script "nohup " & quoted form of scbrowserPath & " serve --port " & serverPort & " > " & quoted form of logPath & " 2>&1 & echo $! > " & quoted form of pidPath
	repeat 180 times
		if my serverUp() then return true
		delay 1
	end repeat
	return false
end startServer

on stopServer()
	set cmd to "if [ -f " & quoted form of pidPath & " ]; then kill $(cat " & quoted form of pidPath & ") 2>/dev/null; rm -f " & quoted form of pidPath & "; fi; "
	set cmd to cmd & "pids=$(/usr/sbin/lsof -ti tcp:" & serverPort & " -sTCP:LISTEN); if [ -n \"$pids\" ]; then kill $pids; fi; true"
	do shell script cmd
end stopServer

on run
	if my serverUp() then
		set answer to button returned of (display dialog ("The gene browser is running at" & return & (my appURL())) buttons {"Stop", "Open"} default button "Open" with title "Gene Browser")
		if answer is "Stop" then
			my stopServer()
			display dialog "The gene browser has stopped." buttons {"OK"} default button "OK" giving up after 3 with title "Gene Browser"
			return
		end if
	else
		set progress total steps to -1
		set progress description to "Starting the gene browser..."
		set progress additional description to "It opens in your browser when ready (about 10-30 seconds)."
		if not (my startServer()) then
			display dialog ("The gene browser did not start. Details are in:" & return & logPath) buttons {"OK"} default button "OK" with icon stop with title "Gene Browser"
			return
		end if
	end if
	open location (my appURL())
end run
