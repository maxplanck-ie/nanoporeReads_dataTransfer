#!/usr/bin/env python
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import smtplib

#send an email:
def sendEmail(body, subject, fromEmail, toEmail, emailHost):
    # set up the SMTP server
    s = smtplib.SMTP(host=emailHost)
    msg = MIMEMultipart() # create a message

    email = MIMEText(body, 'plain')

    # setup the parameters of the message
    msg['From']=fromEmail
    msg['To']=toEmail
    msg['Subject']=subject

    # add in the message body
    msg.attach(email)

    # send the message via the server set up earlier.
    s.send_message(msg)

    del msg
