#!/usr/bin/env python3

__author__ = "Claude 4.0 Sonnet"
__date__ = "June 2, 2025"

"""
Simple script to open the interactive debris cloud HTML visualization in browser.
"""

import webbrowser
import os
import sys

def open_interactive_plot():
    """Open the interactive_cloud.html file in the default browser."""
    
    # Get the path to the HTML file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    html_file = os.path.join(script_dir, "interactive_cloud.html")
    
    # Check if file exists
    if not os.path.exists(html_file):
        print("Error: interactive_cloud.html not found!")
        print("Please run the main simulation first to generate the visualization.")
        return False
    
    # Convert to absolute path
    html_file = os.path.abspath(html_file)
    
    try:
        # Open in default browser
        webbrowser.open(f"file://{html_file}")
        print(f"Opening interactive visualization: {html_file}")
        return True
    except Exception as e:
        print(f"Error opening browser: {e}")
        return False

if __name__ == "__main__":
    open_interactive_plot()
