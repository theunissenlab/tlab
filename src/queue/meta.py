"""SQLAlchemy Metadata and Session object"""
from sqlalchemy import MetaData, create_engine, __version__ as version
from sqlalchemy.orm import scoped_session, sessionmaker

__all__ = ['engine', 'metadata', 'Session']
    
engine = create_engine('mysql://david:1997nine@sql.jlg.berkeley.edu:3306')

session_keywords = {'autoflush':True,
                    'bind':engine,
                    'autocommit':False}

Session = scoped_session(sessionmaker(**session_keywords))

metadata = MetaData(engine)