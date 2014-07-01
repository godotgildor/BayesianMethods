import os
import subprocess
import argparse
import json

import pytumblr

def get_arg_parser():
    ap = argparse.ArgumentParser(description='Converts an ipython notebook into html and then posts to Tumblr.')
    ap.add_argument('-n', '--ipython_notebook',
                    help='The ipython notebook file', required=True)

    ap.add_argument('-o', '--oauth_credentials',
                    help='Json file with oauth credentials.', required=True)

    ap.add_argument('-b', '--blog_url',
                    help='Tumblr blog URL.', required=True)

    ap.add_argument('-t', '--blog_title',
                    help='Tumblr blog title.', required=True)

    return ap.parse_args()

if __name__ == '__main__':
    args = get_arg_parser()

    # First convert the ipython notebook to a webpage
    if os.path.splitext(args.ipython_notebook)[-1] != '.html':
        cmd = ['ipython', 'nbconvert', '--template', 'full', args.ipython_notebook]
        subprocess.check_call(cmd)
        ipython_notebook_webpage = os.path.splitext(args.ipython_notebook)[0] + '.html'
    else:
        ipython_notebook_webpage = args.ipython_notebook

    # Now let's post to tumblr
    with open(args.oauth_credentials) as fh:
        tumblr_credentials = json.loads(fh.read())

    tumblr_client = pytumblr.TumblrRestClient(consumer_key=tumblr_credentials['consumer_key'],
                                              consumer_secret=tumblr_credentials['consumer_secret'],
                                              oauth_token=tumblr_credentials['oauth_token'],
                                              oauth_secret=tumblr_credentials['oauth_secret'])

    with open(ipython_notebook_webpage) as fh:
        html_text = fh.read()
        html_text = html_text.replace('<h1', '<h2')
        html_text = html_text.replace('/h1>', '/h2>')
        response =  tumblr_client.create_text(args.blog_url, body=html_text,
                                              title=args.blog_title, format='html')

        if 'id' in response:
            print 'Blog updated. Post Id: {0}'.format(response['id'])
        else:
            print 'Error posting to blog!'
            print 'Status: {0}'.format(response['meta']['status'])
            print 'Message: {0}'.format(response['meta']['msg'])
            sys.exit(1)
